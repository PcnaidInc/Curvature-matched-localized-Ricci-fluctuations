"""Re-evaluate R13C's normalized-transfer bound with a NEW enclosed preparation C."""
import sys,json,time,hashlib
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'vendor'))
import conditional_budget_R13C as transfer

def run(preparation:Path,active:Path,out:Path,nchi=16384,dps=65):
    if out.exists():raise FileExistsError(out)
    p=json.loads(preparation.read_text())
    if not all(p['checks'].values()) or p['old_R13_C_is_input']:raise ValueError('Missing new preparation evidence')
    transfer.C_REPORTED=p['C_integer']  # Explicit parameterization of copied, unmodified algebra.
    result=transfer.run(active,nchi,dps)
    result['run']='R93R14'
    result['premise']={'C':p['C_integer'],'source':str(preparation.name),
        'preparation_summary_sha256':hashlib.sha256(preparation.read_bytes()).hexdigest(),
        'preparation_enclosure_newly_executed':True,'old_R13_C_is_assumed':False,
        'mathematical_status':'analytic bound with interval-evaluated constants for the explicitly defined smooth SPV-1 equation; software is not formally verified'}
    result['unbounded_terms']=['comparison-width retarded response minus R12 local-order9 trial',
                              'inherited numerical geometry, integration, finite-sum global errors']
    result['certifies_historical_field_arrays']=False
    out.write_text(json.dumps(result,indent=2))
    print(result['conditional_component_bound_upper']);return result
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--preparation',type=Path,required=True);ap.add_argument('--active',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--dps',type=int,default=65);a=ap.parse_args()
    run(a.preparation,a.active,a.out,dps=a.dps)
