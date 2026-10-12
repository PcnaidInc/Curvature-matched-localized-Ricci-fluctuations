"""Only expose the bounded component, never silently a physical continuum value."""
import json
from pathlib import Path

def component(path:Path, *, acknowledge_smooth_model:bool=False):
    if not acknowledge_smooth_model:
        raise ValueError('Acknowledge the stated smooth SPV-1 model and interval-software assumptions')
    d=json.loads(Path(path).read_text())
    if not d.get('preparation_evidence_complete'):
        raise ValueError('New preparation evidence has not passed all declared checks')
    return {'scope':'state-width replacement only', 'upper':d['component_upper_outward_decimal'],
            'model':d['model'], 'remaining_errors':d['unbounded_errors']}

def physical_scalar(*args,**kwargs):
    raise RuntimeError('Unavailable: reference-gap and historical numerical/global errors remain unbounded')
