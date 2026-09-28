from pathlib import Path
import json,numpy as np
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1]
V=json.loads((R/'results/FINAL_VALIDATION.json').read_text());assert V['all_primary_pass']
d=np.genfromtxt(R/'results/COUPLING_COMPARISON.csv',delimiter=',',names=True,dtype=None,encoding='utf8')
sel=d[np.isclose(d['center'],.6)&np.isclose(d['delta_rad'],.25)]
fig,ax=plt.subplots(figsize=(6.2,4.2))
ax.plot(sel['xi'],100*sel['trace_rms_vs_minimal'],marker='o',label='Trace / Ricci-scalar spread')
ax.plot(sel['xi'],100*sel['observer_rms_vs_minimal'],marker='s',label='Observer-direction nontrace spread')
ax.set_xticks([0,1/12,1/6,1/4],['0\nMinimal','1/12','1/6\nConformal','1/4'])
ax.set_xlabel('Scalar curvature coupling (same state and sampling)')
ax.set_ylabel('Standard deviation / minimal value (%)')
ax.set_title('Zero trace does not mean zero nontrace fluctuation')
ax.legend(fontsize=9);fig.tight_layout()
fig.savefig(R/'figures/coupling_comparison.pdf');fig.savefig(R/'figures/coupling_comparison.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(6.2,4.2))
for xi,label in [(0,'Minimal coupling'),(1/6,'Conformal coupling')]:
 p=d[np.isclose(d['xi'],xi)&(d['delta_rad']>0)&(d['delta_rad']<.3)]
 ax.plot(p['center'],p['observer_percent_at_epsilon1e9'],marker='o',label=label)
ax.set_xlim(.62,.28);ax.set_yscale('log');ax.set_xlabel('r / M — farther inward to the right')
ax.set_ylabel('Observer-direction Ricci spread /\nclassical curvature scale (%)')
ax.set_title('The nontrace contribution remains nonzero\nIllustrative gravitational strength ε = 10⁻⁹')
ax.legend(fontsize=9);fig.tight_layout();fig.savefig(R/'figures/nontrace_inward.pdf');fig.savefig(R/'figures/nontrace_inward.png',dpi=180);plt.close(fig)
print('Saved two figures from the completed, validated coupling-control tables.')
