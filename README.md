# Human activity analysis & PT326 identification and control (MATLAB/Simulink)

Two laboratory projects for **Algoritmos de Diagnóstico e de Auto-Regulação** (Diagnosis and Self-Regulation Algorithms), University of Coimbra, Erasmus 2019/2020.
Team: Ferdinando Catania, Giulio Steyde, Hugo Sebastião.

## 1 · Human activity analysis from smartphone sensors — `01-activity-analysis/`

Tri-axial accelerometer and gyroscope signals (50 Hz, smartphone on the waist) from volunteers performing **12 activities**: dynamic (walking, upstairs, downstairs), static (sitting, standing, lying) and transitions.

Pipeline in `assignment_1.m`:
- **ICA** to separate sources, **FIR low-pass** filtering (`firlow2.mat`) and **discrete wavelet decomposition**;
- **feature extraction** in time and frequency domain, statistical selection with **Kruskal–Wallis**;
- dimensionality reduction with **PCA** and **MDS**;
- unsupervised **clustering** (k-medoids, hierarchical), evaluated with silhouette coefficient and external indices against the true labels.

Report: [`docs/report.pdf`](01-activity-analysis/docs/report.pdf). The raw dataset files (~30 MB) are not included.

## 2 · System identification and digital control of the Feedback PT326 — `02-pt326-identification-control/`

The PT326 process trainer heats an air flow; input = heater voltage, output = outlet air temperature (thermistor).

- **Pre-processing**: removal of 1st/2nd-order trends, conversion volt ↔ °C (`data/`).
- **Identification**: ARX models swept over poles (0–20), zeros (0–20) and delay (0–10), selected with **AIC** and **MDL** criteria.
- **Control**: design of a **three-block controller** and a **deadbeat controller**, both continuous and discrete, simulated in Simulink (`*.slx`) for reference tracking and disturbance rejection.

![ARX model selection](02-pt326-identification-control/figures/aic%20criteria.png)

Report: [`docs/report.pdf`](02-pt326-identification-control/docs/report.pdf).

## Skills

MATLAB · Simulink · biosignal processing · ICA / PCA / MDS · wavelets · clustering · system identification (ARX) · digital control (deadbeat)
