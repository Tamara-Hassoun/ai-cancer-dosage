from expert_engine import ColonCancerExpert, Patient
from collections import defaultdict


def get_input(prompt, type_=str):
    while True:
        try:
            value = input(prompt)
            return type_(value) if value != "" else None
        except ValueError:
            print("❌ Invalid input. Please enter a valid", type_.__name__)


def main():
    print("\n🧠 Colorectal Cancer Expert System – Mutation-Based Therapy Assistant\n")

    patient_id = get_input("patient_id (numeric ID): ", float)

    kras_status = get_input("KRAS_status (Mild/Moderate/Severe): ", str)
    kras_impact = get_input("KRAS_impact (LOW/MODERATE/HIGH): ", str)
    kras_change = get_input("KRAS_protein_change (e.g. p.G12D): ", str)

    braf_status = get_input("BRAF_status (Mild/Moderate/Severe): ", str)
    braf_impact = get_input("BRAF_impact (LOW/MODERATE/HIGH): ", str)
    braf_change = get_input("BRAF_protein_change (e.g. p.V600E): ", str)

    msi = get_input("MSI_status (numeric): ", float)
    tmb = get_input("total_perMB (Tumor Mutation Burden): ", float)

    treatment_type = get_input(
        "treatments_type (e.g. Pharmaceutical Therapy, NOS): ", str
    )
    num_cycles = get_input("NUMBER_OF_CYCLES: ", float)
    therapy_ongoing = get_input("THERAPY_ONGOING (Yes/No): ", str)

    stage = get_input("diagnosis_stage (e.g. Stage IIA): ", str)

    vital_status = get_input("vital_status (Alive/Dead): ", str)
    survival_months = get_input("survival_months (numeric): ", float)

    print("\n📊 Processing patient data...\n")

    engine = ColonCancerExpert()
    engine.reset()
    engine.declare(
        Patient(
            KRAS_status=kras_status,
            KRAS_impact=kras_impact,
            KRAS_protein_change=kras_change,
            BRAF_status=braf_status,
            BRAF_impact=braf_impact,
            BRAF_protein_change=braf_change,
            MSI_status=msi,
            total_perMB=tmb,
            treatments_type=treatment_type,
            NUMBER_OF_CYCLES=num_cycles,
            THERAPY_ONGOING=therapy_ongoing,
            diagnosis_stage=stage,
            vital_status=vital_status,
            survival_months=survival_months,
            patient_id=patient_id,
        )
    )
    engine.run()
    engine.print_final_report()


if __name__ == "__main__":
    main()
