"""Bounded HTTP requests, immutable response evidence and independent source statuses."""
import json,time,traceback
from pathlib import Path
import requests
from .provenance import now,sha256,write_json

class RetrievalError(RuntimeError):
    def __init__(self,message,attempts=0,http_status=None):
        super().__init__(message);self.attempts=attempts;self.http_status=http_status


def fetch(url,folder,identifiers=None,params=None):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=False)
    start=now();last=None;status=None
    for attempt in range(1,4):
        try:
            r=requests.get(url,params=params,timeout=(30,120))
            status=r.status_code;raw=folder/f'response_{attempt}.bin';raw.write_bytes(r.content)
            headers={k:v for k,v in r.headers.items() if k.lower() not in ('set-cookie','cookie','authorization')}
            write_json(folder/f'headers_{attempt}.json',headers)
            if r.status_code>=400:
                last=f'HTTP {r.status_code}: {r.text[:200]}'
                if r.status_code not in (408,429,500,502,503,504):break
                delay=min(30,float(r.headers.get('Retry-After',2**(attempt-1))))
            else:
                meta={'source':url.split('/')[2],'database_version':r.headers.get('x-uniprot-release'),
                      'release_date':r.headers.get('x-uniprot-release-date'),'retrieved_at':now(),'timezone':'UTC',
                      'endpoint':r.url,'query_method':'GET','identifiers':identifiers or [],'parameters':params or {},
                      'status':'SUCCESS','http_status':status,'attempts':attempt,'raw_response':str(raw.resolve()),
                      'sha256':sha256(raw),'software':{'requests':requests.__version__},'warnings':[],'errors':[]}
                write_json(folder/'source_manifest.json',meta)
                return r,meta
        except requests.RequestException as e:last=str(e);delay=2**(attempt-1)
        if attempt<3:time.sleep(delay)
    error={'status':'FAILED','error_type':'HTTPError' if status else 'NetworkError','http_status':status,
           'message':last,'attempts':attempt,'started_at':start,'failed_at':now(),'affected_identifiers':identifiers or [],'endpoint':url}
    write_json(folder/'error.json',error);write_json(folder/'source_manifest.json',error)
    raise RetrievalError(last,attempt,status)


def independent_sources(names,callbacks,base,identifiers):
    """Every requested callback is attempted even if every earlier source failed."""
    outcomes=[]
    for name in names:
        start=now();directory=Path(base)/name;directory.mkdir(parents=True,exist_ok=True)
        try:
            result=callbacks[name](directory)
            status=result.get('status','SUCCESS')
            if status not in ('SUCCESS','PARTIAL','FAILED','NOT_RUN'):raise ValueError('Invalid evidence status')
            outcome={'source':name,'status':status,'records':result.get('records',0),'reason':result.get('reason',''),
                     'database_version':result.get('database_version',''),'retrieved_at':result.get('retrieved_at',''),
                     'started_at':start,'finished_at':now(),'error_type':'','http_status':'','attempts':result.get('attempts',0)}
        except Exception as e:
            error={'source':name,'status':'FAILED','error_type':type(e).__name__,'http_status':getattr(e,'http_status',None),
                   'message':str(e),'attempts':getattr(e,'attempts',0),'started_at':start,'failed_at':now(),
                   'affected_identifiers':list(identifiers),'traceback':traceback.format_exc()}
            write_json(directory/'error.json',error)
            write_json(directory/'source_manifest.json',error)
            outcome={'source':name,'status':'FAILED','records':0,'reason':str(e),'database_version':'','retrieved_at':'',
                     'started_at':start,'finished_at':now(),'error_type':type(e).__name__,'http_status':getattr(e,'http_status',''),'attempts':getattr(e,'attempts',0)}
        outcomes.append(outcome)
        print(f"[Stage 05][{name}] {outcome['status']}: {outcome['reason'] or str(outcome['records'])+' records'}",flush=True)
        if outcome['status']=='FAILED':print('No usable evidence from this source; continuing independent sources.',flush=True)
    return outcomes
