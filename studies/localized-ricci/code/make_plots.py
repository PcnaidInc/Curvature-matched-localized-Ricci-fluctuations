from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def main():
 d=np.genfromtxt(ROOT/'results/INWARD_LOCALIZED_RESULTS.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
 fig,ax=plt.subplots(figsize=(8,4.8))
 for obs,label in [('trace','Ricci scalar projection'),('time_projection','Observer-direction Ricci projection')]:
  r=d[(d['observable']==obs)&(d['uniform_sphere']==0)&(d['delta_rad']<.3)]
  ax.plot(r['center'],r['percent_example'],marker='o',label=label)
 ax.set_xlim(.62,.28);ax.set_xlabel('r / M — farther inward to the right');ax.set_ylabel('Standard deviation / classical curvature scale (%)');ax.set_title('Localized contributions grow inward, but remain small\nAll chosen measurement widths track the curvature scale; ε = 10⁻⁹');ax.legend();fig.tight_layout();fig.savefig(ROOT/'figures/inward_gravitational_comparison.png',dpi=180);plt.close(fig)
 fig,ax=plt.subplots(figsize=(8,4.8))
 for obs,label in [('trace','Ricci scalar source'),('time_projection','Observer-direction source')]:
  r=d[(d['observable']==obs)&(d['uniform_sphere']==0)&(d['delta_rad']<.3)]
  ax.plot(r['center'],r['source_to_reference_RMS'],marker='o',label=label)
 ax.axhline(1,linestyle='--');ax.set_xlim(.62,.28);ax.set_xlabel('r / M — farther inward to the right');ax.set_ylabel('Calculated source spread / matched frozen-vacuum reference');ax.set_title('The ordinary-vacuum comparison closely tracks both measurements\nThis is not a ratio to classical gravity');ax.legend();fig.tight_layout();fig.savefig(ROOT/'figures/inward_reference_control.png',dpi=180);plt.close(fig)
if __name__=='__main__':main()
