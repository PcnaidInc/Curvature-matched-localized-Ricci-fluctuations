from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/blackholes-optical-mpl')
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];(ROOT/'figures').mkdir(exist_ok=True)
a=np.genfromtxt(ROOT/'results/FLAT_OPTICAL_REFERENCE.csv',delimiter=',',names=True)
a=a[np.isclose(a['b'],4)]
fig,ax=plt.subplots(figsize=(7,4.3));ax.plot(a['xi'],a['rms_over_G_hbar'],marker='o');ax.set_xlabel('Scalar field–curvature coupling ξ');ax.set_ylabel('Fractional beam-radius RMS coefficient');ax.set_title('Flat-vacuum reference only\nThe conformal case retains a nonzero optical response');fig.tight_layout();fig.savefig(ROOT/'figures/flat_optical_reference.png',dpi=180);plt.close(fig)
a=np.genfromtxt(ROOT/'results/optical_kernel.csv',delimiter=',',names=True)
fig,ax=plt.subplots(figsize=(7,4.3));ax.plot(a['affine_lambda'],a['kernel']);ax.set_xlabel('Affine parameter λ / M');ax.set_ylabel('Causal optical-response weight / M');ax.set_title('Tested Schwarzschild response kernel\nA manufactured-source check, not a quantum-noise prediction');fig.tight_layout();fig.savefig(ROOT/'figures/causal_kernel.png',dpi=180);plt.close(fig)
