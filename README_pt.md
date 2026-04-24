# Semi-Derivadas em Voltametria — Implementação Numérica

Cálculo numérico de derivadas fracionárias (ordem 1/2) pelos métodos de **Grünwald–Letnikov** e **Riemann–Liouville**, aplicado ao perfil de concentração de superfície em voltametria reversível.

> Código suplementar ao artigo: *Understanding the Mathematics of Voltammetry: A Step-by-Step Guide* — Chioquetti, Vital & Serrano, IQ-USP (2026).

---

## Dependências

```bash
pip install numpy scipy matplotlib
```

| Pacote | Versão mínima | Uso |
|--------|--------------|-----|
| `numpy` | 1.20 | Operações vetoriais, `np.gradient` |
| `scipy` | 1.7 | Função Gamma (`scipy.special.gamma`) |
| `matplotlib` | 3.4 | Geração dos gráficos |

> `pandas` está importado no script mas não é utilizado — pode ser removido sem impacto.

---

## Como Executar

```bash
python voltammetry_semiderivative.py
```

Uma janela do matplotlib será aberta com dois painéis: o perfil sigmoidal e as semi-derivadas calculadas pelos dois métodos.

---

## Funções

### `gl_fractional_derivative(x, y, alpha=0.5)`

Derivada fracionária de ordem `alpha` pelo método de **Grünwald–Letnikov**.

**Parâmetros**

| Nome | Tipo | Descrição |
|------|------|-----------|
| `x` | array-like | Variável independente. **Deve ser igualmente espaçada.** |
| `y` | array-like | Valores da função a diferenciar. |
| `alpha` | float | Ordem da operação. `alpha=0.5` → semi-derivada. `alpha=-0.5` → semi-integral. Padrão: `0.5`. |

**Retorno:** `np.ndarray` com os valores da derivada fracionária, mesmo comprimento que `x` e `y`.

**Raises:** `ValueError` se `x` não for igualmente espaçado.

**Exemplo de uso**

```python
import numpy as np
x = np.linspace(0, 10, 1000)
y = np.sin(x)
dy_half = gl_fractional_derivative(x, y, alpha=0.5)
```

**Funcionamento interno**

O método aproxima a derivada fracionária pela soma discreta:

$$\frac{d^\alpha y}{dx^\alpha}\bigg|_{x_i} \approx \frac{1}{h^\alpha} \sum_{k=0}^{i} (-1)^k \binom{\alpha}{k} y_{i-k}$$

Os coeficientes binomiais generalizados são calculados recursivamente:

```python
coeffs[0] = 1.0
coeffs[k] = coeffs[k-1] * (alpha - k + 1) / k
```

> **Atenção:** complexidade $O(n^2)$. Para séries com mais de ~50.000 pontos, considere vetorizar a soma interna com `numpy` ou usar FFT.

---

### `rl_fractional_derivative(x, y, alpha=0.5)`

Derivada fracionária de ordem `alpha` pelo método de **Riemann–Liouville**.

**Parâmetros**

| Nome | Tipo | Descrição |
|------|------|-----------|
| `x` | array-like | Variável independente. **Deve ser igualmente espaçada.** |
| `y` | array-like | Valores da função a diferenciar. |
| `alpha` | float | Ordem da operação. Padrão: `0.5`. |

**Retorno:** `np.ndarray` com os valores da derivada fracionária.

**Raises:** `ValueError` se `x` não for igualmente espaçado.

**Exemplo de uso**

```python
dy_half_rl = rl_fractional_derivative(x, y, alpha=0.5)
```

**Funcionamento interno**

Primeiro calcula a integral fracionária (quadratura retangular), depois diferencia com `np.gradient`:

```python
# Integral fracionária
I[j] = (h / gamma(alpha)) * sum(y[k] * (x[j] - x[k])**(alpha-1) for k != j)

# Diferenciação numérica
rl_frac_deriv = np.gradient(I, h)
```

> O termo `k == j` é excluído porque o integrando $(x_j - x_k)^{\alpha-1}$ é singular em $k = j$ para $\alpha < 1$.

---

## Parâmetros Físicos do Script

Os parâmetros a seguir definem o sistema eletroquímico simulado e podem ser ajustados diretamente no script:

| Variável | Valor padrão | Unidade | O que controla |
|----------|-------------|---------|----------------|
| `n` | `1` | — | Número de elétrons transferidos |
| `T` | `298` | K | Temperatura (padrão: 25 °C) |
| `F` | `96485` | C/mol | Constante de Faraday (não alterar) |
| `R` | `8.314` | J/(K·mol) | Constante dos gases (não alterar) |
| `scan_rate` | `0.05` | V/s | Velocidade de varredura do potencial |
| `delta_E` | `0.5` | V | $E_i - E^\circ$: afasta o potencial inicial do $E^\circ$ |

As variáveis adimensionais `a` e `b` são derivadas automaticamente:

```python
a = (n * F * scan_rate) / (R * T)   # taxa de varredura adimensional [s⁻¹]
b = (n * F * delta_E) / (R * T)     # potencial inicial adimensional
theta = -a * time + b               # potencial adimensional θ(t)
```

A função de entrada para as semi-derivadas é a **sigmoidal de Nernst**:

```python
y = 1 / (1 + np.exp(theta))
```

---

## Saída

O script produz uma figura com dois painéis:

**Painel esquerdo — Perfil de concentração**
- Curva sigmoidal $f(\theta) = 1/(1+e^\theta)$ em função de $\theta$
- Representa a concentração normalizada de Ox na superfície do eletrodo

**Painel direito — Semi-derivadas (perfil de corrente)**
- $-d^{1/2}f / d(at)^{1/2}$ em função de $\theta$
- Duas curvas sobrepostas: GL (tracejada) e RL (tracejada)
- Representa o voltamograma adimensional: pico em $\theta \approx -1.11$, altura $\approx -0.446$

As curvas GL e RL devem coincidir — qualquer divergência expressiva indica problema numérico (espaçamento irregular, `alpha` fora de [0,1], série muito curta).

---

## Adaptando para Seus Dados

Para aplicar as funções a dados próprios, substitua `x` e `y` pelos seus vetores (igualmente espaçados):

```python
import numpy as np

# Seus dados experimentais ou simulados
E = np.linspace(-0.5, 0.5, 5000)   # potencial [V]
concentration_profile = ...         # C(E) normalizada

semi_deriv = gl_fractional_derivative(E, concentration_profile, alpha=0.5)
```

Para calcular a semi-integral (ordem $-1/2$), passe `alpha=-0.5`:

```python
semi_integral = gl_fractional_derivative(E, current_data, alpha=-0.5)
```

---

## Limitações Conhecidas

- **Espaçamento irregular:** ambas as funções lançam `ValueError`. Interpole antes com `np.interp` ou `scipy.interpolate`.
- **Efeito de borda (GL):** os primeiros ~10–20 pontos têm menor precisão por memória fracionária insuficiente.
- **Singularidade (RL):** a exclusão do ponto $k = j$ introduz erro crescente para `alpha` próximo de 0. Para `alpha = 0.5` o erro é aceitável para visualização.
- **`alpha` fora de (0, 1):** o código aceita qualquer valor real, mas resultados para `alpha > 1` ou `alpha < 0` não foram validados neste contexto.

---

## Referências

- Chioquetti, R. A. L.; Vital, J. V. S.; Serrano, S. H. P. *Understanding the Mathematics of Voltammetry: A Step-by-Step Guide.* IQ-USP, 2026.
- Oldham, K. B.; Spanier, J. *Fractional Calculus.* Academic Press, 1974.
- Oldham, K. B. *Fractional differential equations in electrochemistry.* Adv. Eng. Softw., 41, 9–12, 2010.
