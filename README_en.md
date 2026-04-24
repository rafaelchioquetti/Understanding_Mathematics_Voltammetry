# Fractional Semi-Derivatives in Voltammetry — Numerical Implementation

Numerical computation of fractional derivatives (order 1/2) using the **Grünwald–Letnikov** and **Riemann–Liouville** methods, applied to the surface concentration profile of a reversible voltammetric system.

> Supplementary code for: *Understanding the Mathematics of Voltammetry: A Step-by-Step Guide* — Chioquetti, Vital & Serrano, IQ-USP (2026).

---

## Dependencies

```bash
pip install numpy scipy matplotlib
```

| Package | Min. version | Usage |
|---------|-------------|-------|
| `numpy` | 1.20 | Array operations, `np.gradient` |
| `scipy` | 1.7 | Gamma function (`scipy.special.gamma`) |
| `matplotlib` | 3.4 | Plot generation |

> `pandas` is imported in the script but not used — it can be safely removed.

---

## Running the Script

```bash
python voltammetry_semiderivative.py
```

A matplotlib window will open showing two panels: the sigmoidal concentration profile and the semi-derivatives computed by both methods.

---

## Functions

### `gl_fractional_derivative(x, y, alpha=0.5)`

Fractional derivative of order `alpha` using the **Grünwald–Letnikov** method.

**Parameters**

| Name | Type | Description |
|------|------|-------------|
| `x` | array-like | Independent variable. **Must be evenly spaced.** |
| `y` | array-like | Function values to differentiate. |
| `alpha` | float | Order of the operation. `alpha=0.5` → semi-derivative. `alpha=-0.5` → semi-integral. Default: `0.5`. |

**Returns:** `np.ndarray` with fractional derivative values, same length as `x` and `y`.

**Raises:** `ValueError` if `x` is not evenly spaced.

**Usage example**

```python
import numpy as np
x = np.linspace(0, 10, 1000)
y = np.sin(x)
dy_half = gl_fractional_derivative(x, y, alpha=0.5)
```

**How it works**

The method approximates the fractional derivative by the discrete sum:

$$\frac{d^\alpha y}{dx^\alpha}\bigg|_{x_i} \approx \frac{1}{h^\alpha} \sum_{k=0}^{i} (-1)^k \binom{\alpha}{k} y_{i-k}$$

The generalized binomial coefficients are computed recursively:

```python
coeffs[0] = 1.0
coeffs[k] = coeffs[k-1] * (alpha - k + 1) / k
```

> **Performance note:** $O(n^2)$ complexity. For series longer than ~50,000 points, consider vectorizing the inner sum with `numpy` or using an FFT-based approach.

---

### `rl_fractional_derivative(x, y, alpha=0.5)`

Fractional derivative of order `alpha` using the **Riemann–Liouville** method.

**Parameters**

| Name | Type | Description |
|------|------|-------------|
| `x` | array-like | Independent variable. **Must be evenly spaced.** |
| `y` | array-like | Function values to differentiate. |
| `alpha` | float | Order of the operation. Default: `0.5`. |

**Returns:** `np.ndarray` with fractional derivative values.

**Raises:** `ValueError` if `x` is not evenly spaced.

**Usage example**

```python
dy_half_rl = rl_fractional_derivative(x, y, alpha=0.5)
```

**How it works**

First computes the fractional integral (rectangular quadrature), then differentiates with `np.gradient`:

```python
# Fractional integral
I[j] = (h / gamma(alpha)) * sum(y[k] * (x[j] - x[k])**(alpha-1) for k != j)

# Numerical differentiation
rl_frac_deriv = np.gradient(I, h)
```

> The `k == j` term is excluded because the integrand $(x_j - x_k)^{\alpha-1}$ is singular at $k = j$ for $\alpha < 1$.

---

## Physical Parameters

The following parameters define the simulated electrochemical system and can be freely adjusted in the script:

| Variable | Default | Unit | Controls |
|----------|---------|------|----------|
| `n` | `1` | — | Number of electrons transferred |
| `T` | `298` | K | Temperature (default: 25 °C) |
| `F` | `96485` | C/mol | Faraday constant (do not change) |
| `R` | `8.314` | J/(K·mol) | Gas constant (do not change) |
| `scan_rate` | `0.05` | V/s | Potential scan rate |
| `delta_E` | `0.5` | V | $E_i - E^\circ$: sets the initial potential relative to $E^\circ$ |

The dimensionless variables `a` and `b` are derived automatically:

```python
a = (n * F * scan_rate) / (R * T)   # dimensionless scan rate [s⁻¹]
b = (n * F * delta_E) / (R * T)     # dimensionless initial potential
theta = -a * time + b               # dimensionless potential θ(t)
```

The input function for the semi-derivatives is the **Nernst sigmoidal**:

```python
y = 1 / (1 + np.exp(theta))
```

---

## Output

The script produces a figure with two panels:

**Left panel — Concentration profile**
- Sigmoidal curve $f(\theta) = 1/(1+e^\theta)$ as a function of $\theta$
- Represents the normalized surface concentration of Ox at the electrode

**Right panel — Semi-derivatives (current profile)**
- $-d^{1/2}f / d(at)^{1/2}$ as a function of $\theta$
- Two overlapping dashed curves: GL and RL
- Represents the dimensionless voltammogram: peak at $\theta \approx -1.11$, height $\approx -0.446$

The GL and RL curves should overlap — any significant divergence points to a numerical issue (uneven spacing, `alpha` outside [0,1], or too short a series).

---

## Adapting to Your Own Data

To apply the functions to your own data, replace `x` and `y` with your (evenly spaced) arrays:

```python
import numpy as np

# Your experimental or simulated data
E = np.linspace(-0.5, 0.5, 5000)   # potential [V]
concentration_profile = ...         # normalized C(E)

semi_deriv = gl_fractional_derivative(E, concentration_profile, alpha=0.5)
```

To compute the semi-integral (order $-1/2$), pass `alpha=-0.5`:

```python
semi_integral = gl_fractional_derivative(E, current_data, alpha=-0.5)
```

---

## Known Limitations

- **Uneven spacing:** both functions raise `ValueError`. Interpolate first using `np.interp` or `scipy.interpolate`.
- **Edge effect (GL):** the first ~10–20 points have lower accuracy due to insufficient fractional memory.
- **Singularity (RL):** excluding the `k == j` point introduces growing error as `alpha` approaches 0. For `alpha = 0.5`, the error is acceptable for visualization purposes.
- **`alpha` outside (0, 1):** the code accepts any real value, but results for `alpha > 1` or `alpha < 0` have not been validated in this context.

---

## References

- Chioquetti, R. A. L.; Vital, J. V. S.; Serrano, S. H. P. *Understanding the Mathematics of Voltammetry: A Step-by-Step Guide.* IQ-USP, 2026.
- Oldham, K. B.; Spanier, J. *Fractional Calculus.* Academic Press, 1974.
- Oldham, K. B. *Fractional differential equations in electrochemistry.* Adv. Eng. Softw., 41, 9–12, 2010.
