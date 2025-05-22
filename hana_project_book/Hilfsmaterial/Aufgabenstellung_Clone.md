# 1 Aufgabestellung

Das elektromagnetische Feld kann mit Hilfe der Maxwell Gleichungen beschrieben werden.
Für die niederfrequente Anwendung wird die harmonische Wirbelstromgleichung

$$ \sigma j\omega \mathbf{A} + \text{curl}\mu^{-1}\text{curl}\mathbf{A} + \sigma\nabla\varphi = \mathbf{j}_{\text{ext}} \quad \text{für } x \in \Omega_{\sigma>0} \quad \text{(1a)} $$

$$ \text{div}(\sigma j\omega \mathbf{A} + \sigma\nabla\varphi) = \text{div}(\mathbf{j}_{\text{ext}}) \quad \text{für } x \in \Omega_{\sigma>0} \quad \text{(1b)} $$

$$ \text{curl}\mu^{-1}\text{curl}\mathbf{A} = \mathbf{j}_{\text{ext}} \quad \text{für } x \in \Omega_{\sigma=0} \quad \text{(1c)} $$

benutzt, wobei mit $\mathbf{A}$ das Vektorpotential und mit $\varphi$ das skalare Potential bezeichnet werden. Es gilt

$$ \mathbf{B} = \text{curl}\mathbf{A} $$

und

$$ \mathbf{E} = -\sigma j\omega \mathbf{A} - \nabla\varphi. $$

Die zeitabhängige Lösung ist gegeben durch

$$ \mathbf{A}(t, x) = \mathbf{A}(x) \cdot e^{j\omega t} $$

mit $\omega = 2\pi f$ und $j = \sqrt{-1}$.

Die Stromdichte $\mathbf{J}$ kann in zwei Komponenten aufgeteilt werden. Es gilt

$$ \mathbf{J} = \mathbf{J}_{\text{int}} + \mathbf{J}_{\text{ext}}. $$

mit der induzierten Stromdichte $\mathbf{J}_{\text{int}}$ (Wirbelströme) und der externen Stromdichte $\mathbf{J}_{\text{ext}}$.
Die externe Stromdichte ist zum Beispiel gegeben durch eine Stromquelle, mit welcher die Spulen getrieben werden.
Für dreidimensionale Geometrien ist es in der Regel nicht möglich sämtliche Windungen im Mesh aufzulösen. Unter Vernachlässigung von Proximity- und Skineffekt in den Windungen folgt mit $0 < \kappa < 1$ das regularisierte Wirbelstromproblem

$$ \sigma j\omega \mathbf{A} + \text{curl}\mu^{-1}\text{curl}\mathbf{A} + \kappa\mathbf{A} = \mathbf{j}_{\text{ext}} \quad \text{für } x \in \Omega. \quad \text{(2)} $$

Für den Fall, dass eine **Stromquelle** gegeben ist, gilt

$$ \mathbf{J}_{\text{ext}} = \frac{N_c i_c}{A_c} \mathbf{w}(x), $$

mit Anzahl Windungen $N_c$, Strom in den Windungen $i_c$ und der Querschnittsfläche $A_c$.

Für den Fall, dass eine **Spannungsquelle** gegeben ist, folgt aus dem Induktionsgesetz die Kopplung zwischen Strom und Spannung in den Windungen. Es gilt

$$ j\omega\Psi + R_c i_c = u_c. $$

Damit folgt das gekoppelte Problem

$$ \sigma j\omega \mathbf{A} + \text{curl}\mu^{-1}\text{curl}\mathbf{A} + \kappa\mathbf{A} - \frac{N_c}{A_c}\mathbf{w}(x) i_c = 0 \quad \text{für } x \in \Omega. $$

$$ j\omega\Psi + R_c i_c = u_c. \quad \text{(3)} $$

## 1.1 Rotationssymmetrische Gleichung

Für rotationssymmetrische Geometrien können in der Regel die Windungen im Mesh aufgelöst werden (vgl. Abbildung 1). Wir betrachten daher das Problem (1) in Zylinderkoordinaten $(z, r, \phi)$. In diesem Fall gilt

$$ \mathbf{A} = \begin{pmatrix} A_z \\ A_r \\ A_\phi \end{pmatrix} $$

wobei $A_z = A_r = 0$. Das Vektorpotential ist daher gegeben durch die skalare Komponente $A_\phi$. Damit folgt für (1)

$$ \sigma j\omega A_\phi - \left( \frac{\partial}{\partial r} \left( \frac{\mu^{-1}}{r} \frac{\partial}{\partial r} (r A_\phi) \right) + \frac{\partial}{\partial z} \left( \mu^{-1} \frac{\partial A_\phi}{\partial z} \right) \right) = j_\phi. \quad \text{(4)} $$

Für die magnetische Flussdichte gilt in dem Fall

$$ B_r = -\frac{\partial A_\phi}{\partial z} $$

$$ B_z = \frac{1}{r} \frac{\partial (r A_\phi)}{\partial r} \quad \text{(5)} $$

Das elektrische Feld ist orthogonal zum Magnetfeld und zeigt in $\phi$ Richtung. Damit folgt

$$ E_z = -j\omega A_\phi. \quad \text{(6)} $$

Die natürliche Symmetriebedingung auf der Rotationsachse ist gegeben durch

$$ A_\phi = 0. $$

<!-- Abbildung 1: Geometrie -->
*(Description: Abbildung 1 zeigt eine Mesh-Geometrie mit einem roten Leiter und grünen Windungen.)*

Im rotationssymmetrischen Modell, mit den einzelnen Windungen modelliert, folgt für die Stromdichte $j_\phi$

$$ j_\phi = \frac{i_c}{r^2\pi}. $$

Wir setzen $r = 0.001\text{m}$ und $i_c = 10\text{A}$. Die magnetische Permeabilität $\mu$ ist gegeben durch $\mu = \mu_r \cdot \mu_0$. Die relative Permeabilität ist ein Mass für die Feldverstärkungim Material (vgl. Tabelle 1). Vom Gebiet abhängige Konstanten können mit Hilfe von `CoefficientFunction` definiert werden (vgl. Listing 7).

| Material                     | Luft                               | Windungen                          | Kern                               |
|------------------------------|------------------------------------|------------------------------------|------------------------------------|
| el. Leitfähigkeit $\sigma$   | $0 \text{ S/m}$                  | $56 \cdot 10^6 \text{ S/m}$      | $56 \cdot 10^6 \text{ S/m}$      |
| relative mag. Permeabilität $\mu_r$ | $1$                              | $1$                              | $1$                              |
| Wärmeleitfähigkeit $\lambda$ | $0.0262 \text{ W/(m K)}$         | $400 \text{ W/(m K)}$            | $400 \text{ W/(m K)}$            |
Tabelle 1: Materialparameter.

|                              |             |
|------------------------------|-------------|
| Frequenz $f$               | $50 \text{ Hz}$ |
| magnetische Feldkonstante $\mu_0$ | $4\pi 10^{-7}$ |
Tabelle 2: Konstanten

## 1.2 Induktive Wärmequelle

Die induktive Wärmequelle ist gegeben durch die mittleren Wärmeverluste des elektromagnetischen Feldes

$$ q = \frac{1}{2} |\mathbf{J} \cdot \mathbf{E}|= \frac{\sigma}{2} |\mathbf{E}|^2, \quad \text{(7)} $$

wobei mit $\mathbf{J} = \mathbf{J}_{\text{int}} + \mathbf{J}_{\text{ext}}$ die totale Stromdichte bezeichnet sei.

## 1.3 Thermisches Feld

Die parabolische Wärmeleitungsgleichung ist gegeben durch

$$ C_p \rho \dot{T} - \text{div}(\lambda\nabla T) + \nabla T \cdot \mathbf{v} = q, \quad \text{(8)} $$

wobei in der Form auch der konvektive Wärmetransport gegeben durch ein Strömungsfeld $\mathbf{v}$ berücksichtigt wird. Im Folgenden wird dieser Term für das erste vernachlässigt.
Die stationäre Lösung der Wärmeleitungsgleichung (8) genügt der elliptischen Differentialgleichung

$$ -\text{div}(\lambda\nabla T) = q. \quad \text{(9)} $$

**Randbedingungen:**

1.  konstante Temperatur am Rand (Dirichlet Randbedingung)

    $$ T(x) = 0 \quad x \in \partial\Omega $$
2.  perfekte Temperatur Isolation am Rand (Neumann Randbedingung)

    $$ \nabla T \cdot \mathbf{n} = 0 \quad x \in \partial\Omega $$

