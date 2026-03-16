# 🧠 Mutation-Based Therapy Expert System for Colorectal Cancer

## 🎯 Project Goal
An interactive expert system designed to predict **patient response to cancer treatment** based on **molecular, clinical, and therapeutic** features.

The system simulates a **human clinical decision process** using **if-then rules** and knowledge extracted from real patient data.

---

## 🧬 Features Used

- **Mutational Data**
  - KRAS/BRAF status, impact level, protein changes
- **Immune Markers**
  - MSI (Microsatellite Instability)
  - Tumor Mutation Burden (TMB)
- **Therapy History**
  - Type of treatment
  - Number of cycles
  - Ongoing status
- **Clinical Stage & Vital Data**
  - AJCC stage, survival months, vital status
  - Patient ID (used in ML-derived logic)

---

## ⚙️ How to Run

```bash
pip install experta
python main.py
