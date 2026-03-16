from experta import *


class Patient(Fact):
    """Fact: Full clinical and molecular profile for colorectal cancer decision support"""

    pass


class Result(Fact):
    type = Field(str)  
    message = Field(str)  
    certainty_factor = Field(float, default=1.0)  





class ColonCancerExpert(KnowledgeEngine):
    def combine_cf(self, cf1, cf2):
        if cf1 > 0 and cf2 > 0:
            return round(cf1 + cf2 - cf1 * cf2, 2)
        elif cf1 < 0 and cf2 < 0:
            return round(cf1 + cf2 + cf1 * cf2, 2)
        else:
            denominator = 1 - min(abs(cf1), abs(cf2))
            return round((cf1 + cf2) / denominator, 2) if denominator != 0 else 0.0

    def __init__(self):
        super().__init__()
        self.outputs = {
            "Genetic": [],
            "Immune": [],
            "Therapy": [],
            "Survival": [],
            "ML": [],
            "Exceptional": [],
            "Conflict": [],
            "Recommendation": [],
            "Fallback": [],
        }

    @DefFacts()
    def _initial_action(self):
        yield Fact(action="analyze_patient")

    # combine rule

    @Rule(
        AS.f1 << Result(type=MATCH.t, message=MATCH.m, certainty_factor=MATCH.cf1),
        AS.f2 << Result(type=MATCH.t, message=MATCH.m, certainty_factor=MATCH.cf2),
        TEST(lambda f1, f2: f1 != f2),
    )
    def combine_duplicate_results(self, f1, f2, cf1, cf2):
        self.retract(f1)
        combined_cf = self.combine_cf(cf1, cf2)
        self.modify(f2, certainty_factor=combined_cf)

    # ─── conflict cases ───
    @Rule(
        Fact(action="analyze_patient"),
        Patient(MSI_status=P(lambda x: x < 0.5)),
        Patient(msi_cf=P(lambda x: x >= 0.85)),
        Patient(survival_months=P(lambda x: x > 100)),
    )
    def rule_msi_low_but_survives_conflict(self):
        cf_rule = 0.60
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: MSI reported as low with high input confidence, but survival is long — review MSI methodology or patient-specific immunity.",
                certainty_factor=cf_rule,
                source="PMID: 29887279",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(MSI_status=P(lambda x: x < 0.5)),
        Patient(total_perMB=P(lambda x: x < 5)),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 3)),
        Patient(survival_months=P(lambda x: x >= 80)),
    )
    def rule_conflict_nonimmune_but_survives(self):
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Low immune activity and short therapy but survival >6 years — immune-independent tumor suppression likely.",
                certainty_factor=0.65,
                source="PMID: 29887279",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Severe"),
        Patient(BRAF_status="Severe"),
        Patient(MSI_status=P(lambda x: x < 0.5)),
        Patient(total_perMB=P(lambda x: x < 5)),
    )
    def rule_dual_oncogene_low_immune(self):
        cf_rule = 0.70
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Severe KRAS/BRAF with low MSI/TMB — likely immune evasion and high resistance.",
                certainty_factor=cf_rule,
                source="PMID: 37059545",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Severe"),
        Patient(MSI_status=P(lambda x: x < 0.5)),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 3)),
    )
    def rule_combined_kras_lowimmune_shorttherapy(self):
        cf_rule = 0.80
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Severe KRAS + low immune score + short therapy — high risk of recurrence or progression.",
                certainty_factor=cf_rule,
                source="PMID: 35788316",
            )
        )

    # ─── GENETIC GENERAL RULES ───

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_status="Severe", KRAS_impact="HIGH", kras_cf=P(lambda x: x is not None)
        ),
    )
    def rule_kras_severe_high(self, kras_cf):
        cf_rule = 0.95
        cf_conclusion = cf_rule * kras_cf
        self.declare(
            Result(
                type="Genetic",
                message=" No Response: Severe KRAS mutation with high impact is associated with drug resistance.",
                certainty_factor=cf_conclusion,
                source="OncoKB",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            BRAF_status="Severe", BRAF_impact="HIGH", braf_cf=P(lambda x: x is not None)
        ),
    )
    def rule_braf_high_impact(self, braf_cf):
        cf_rule = 0.75
        cf_conclusion = cf_rule * braf_cf
        self.declare(
            Result(
                type="Genetic",
                message=" Stable: BRAF high impact mutation may cause variable response depending on treatment strategy.",
                certainty_factor=cf_conclusion,
                source="NCCN",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Consider combination of BRAF and MEK inhibitors for BRAF high-impact cases.",
                certainty_factor=0.85,
                source="PMID: 26559394",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_status="Mild", KRAS_impact="LOW", kras_cf=P(lambda x: x is not None)
        ),
    )
    def rule_kras_mild(self, kras_cf):
        cf_rule = 0.90
        cf_conclusion = cf_rule * kras_cf
        self.declare(
            Result(
                type="Genetic",
                message=" Response: Mild KRAS mutation with low impact may indicate better drug sensitivity.",
                certainty_factor=cf_conclusion,
                source="OncoKB",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(BRAF_protein_change="p.V600E", braf_cf=P(lambda x: x is not None)),
    )
    def rule_braf_v600e(self, braf_cf):
        cf_rule = 0.95
        cf_conclusion = cf_rule * braf_cf
        self.declare(
            Result(
                type="Genetic",
                message=" Response: BRAF V600E is targetable by BRAF inhibitors like vemurafenib.",
                certainty_factor=cf_conclusion,
                source="PMID: 26559394",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggestion: Use BRAF-targeted therapy such as vemurafenib or encorafenib.",
                certainty_factor=0.95,
                source="NCCN/OncoKB",
            )
        )

    # ─── SPECIFIC PROTEIN MUTATIONS ───

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_protein_change="p.G12V",
            KRAS_status="Severe",
            kras_cf=P(lambda x: x is not None),
        ),
    )
    def rule_kras_g12v(self, kras_cf):
        cf_rule = 0.90
        cf_conclusion = cf_rule * kras_cf
        self.declare(
            Result(
                type="Genetic",
                message=" No Response: KRAS G12V is linked to strong resistance to EGFR-targeted therapies.",
                certainty_factor=cf_conclusion,
                source="OncoKB",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_protein_change="p.Q61L",
            KRAS_impact="MODERATE",
            kras_cf=P(lambda x: x is not None),
        ),
    )
    def rule_kras_q61l(self, kras_cf):
        cf_rule = 0.75
        cf_conclusion = cf_rule * kras_cf
        self.declare(
            Result(
                type="Genetic",
                message=" Stable: KRAS Q61L may lead to partial drug response or stable disease.",
                certainty_factor=cf_conclusion,
                source="GENIE/CBioPortal",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_protein_change="p.G12D", kras_cf=P(lambda x: x is not None)),
    )
    def rule_kras_g12d(self, kras_cf):
        cf_rule = 0.90
        cf_conclusion = cf_rule * kras_cf
        self.declare(
            Result(
                type="Genetic",
                message=" No Response: KRAS G12D mutation generally resists anti-EGFR treatment.",
                certainty_factor=cf_conclusion,
                source="OncoKB",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(BRAF_protein_change="p.K601E", braf_cf=P(lambda x: x is not None)),
        Patient(survival_months=P(lambda x: x >= 24)),
    )
    def rule_braf_k601e_positive(self, braf_cf):
        cf_rule = 0.85
        cf_conclusion = cf_rule * braf_cf
        self.declare(
            Result(
                type="Genetic",
                message=" Response: BRAF K601E with favorable survival suggests potential benefit from BRAF/MEK pathway targeting.",
                certainty_factor=cf_conclusion,
                source="PMID: 34445615",
            )
        )

    # ─── IMMUNE PROFILE & TMB ───

    @Rule(
        Fact(action="analyze_patient"),
        Patient(MSI_status=P(lambda x: x >= 1.0), msi_cf=MATCH.msi_cf),
    )
    def rule_high_msi(self, msi_cf):
        cf_rule = 0.95
        cf_conclusion = cf_rule * msi_cf
        self.declare(
            Result(
                type="Immune",
                message=" Response: High Microsatellite Instability (MSI ≥ 1.0) strongly indicates potential benefit from immunotherapy.",
                certainty_factor=cf_conclusion,
                source="NCCN Guidelines",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Consider initiating checkpoint inhibitors (e.g., pembrolizumab) in MSI-high patients.",
                certainty_factor=0.95,
                source="PMID: 29887279",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            MSI_status=P(lambda x: x < 0.1),
            total_perMB=P(lambda x: x < 2.0),
            msi_cf=MATCH.msi_cf,
        ),
    )
    def rule_low_msi_tmb(self, msi_cf):
        cf_rule = 0.85
        cf_conclusion = cf_rule * msi_cf
        self.declare(
            Result(
                type="Immune",
                message=" No Response: Low MSI and TMB indicate immunologically 'cold' tumor — poor response to immunotherapy.",
                certainty_factor=cf_conclusion,
                source="PMID: 29887279",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(total_perMB=P(lambda x: x > 20), tmb_cf=MATCH.tmb_cf),
    )
    def rule_high_tmb(self, tmb_cf):
        cf_rule = 0.90
        cf_conclusion = cf_rule * tmb_cf
        self.declare(
            Result(
                type="Immune",
                message=" Response: High TMB (>20 mut/Mb) is associated with increased neoantigen burden and enhanced immune recognition.",
                certainty_factor=cf_conclusion,
                source="FDA Pembrolizumab Approval",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Consider immunotherapy (e.g., anti-PD-1) for patients with high TMB.",
                certainty_factor=0.90,
                source="FDA, NCCN",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            MSI_status=P(lambda x: x >= 1.0),
            total_perMB=P(lambda x: x >= 10),
            msi_cf=MATCH.msi_cf,
            tmb_cf=MATCH.tmb_cf,
        ),
    )
    def rule_high_msi_high_tmb(self, msi_cf, tmb_cf):
        cf_rule = 0.95
        cf_conclusion = cf_rule * msi_cf * tmb_cf
        self.declare(
            Result(
                type="Immune",
                message=" Response: Combined MSI-high and TMB-high profile strongly predicts durable immunotherapy response.",
                certainty_factor=cf_conclusion,
                source="NCCN, OncoKB",
            )
        )

    # ─── THERAPY-RELATED RULES ───

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            treatments_type="Pharmaceutical Therapy, NOS",
            NUMBER_OF_CYCLES=P(lambda x: x < 6),
        ),
    )
    def rule_short_pharma_treatment(self):
        cf_rule = 0.80
        self.declare(
            Result(
                type="Therapy",
                message=" No Response: Fewer than 6 cycles of pharmaceutical therapy often lead to suboptimal outcomes.",
                certainty_factor=cf_rule,
                source="Clinical Practice",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Consider optimizing treatment duration or switching to multimodal regimen in poor responders.",
                certainty_factor=0.85,
                source="Clinical Practice",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            treatments_type="Radiation Therapy, NOS",
            NUMBER_OF_CYCLES=P(lambda x: x >= 10),
        ),
    )
    def rule_long_radiation(self):
        cf_rule = 0.75
        self.declare(
            Result(
                type="Therapy",
                message=" Stable: Prolonged radiation therapy (≥10 cycles) may achieve disease stabilization in select patients.",
                certainty_factor=cf_rule,
                source="OncoRadiology Review",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(THERAPY_ONGOING="Yes", NUMBER_OF_CYCLES=P(lambda x: x >= 12)),
    )
    def rule_ongoing_treatment_positive(self):
        cf_rule = 0.85
        self.declare(
            Result(
                type="Therapy",
                message=" Response: Ongoing treatment for ≥12 cycles suggests clinical benefit and disease control.",
                certainty_factor=cf_rule,
                source="Clinical Duration Indicator",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(THERAPY_ONGOING="No", NUMBER_OF_CYCLES=P(lambda x: x <= 3)),
    )
    def rule_early_stop(self):
        cf_rule = 0.80
        self.declare(
            Result(
                type="Therapy",
                message=" No Response: Early discontinuation of therapy (≤3 cycles) may reflect resistance, intolerance, or failure.",
                certainty_factor=cf_rule,
                source="Oncology Practice",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(THERAPY_ONGOING="No", NUMBER_OF_CYCLES=P(lambda x: x <= 3)),
    )
    def rule_short_therapy(self):
        cf_rule = 0.80
        self.declare(
            Result(
                type="Therapy",
                message=" No Response: Early therapy discontinuation (≤3 cycles) is often associated with poor outcomes.",
                certainty_factor=cf_rule,
                source="PMID: 33140315",
            )
        )

    # ─── STAGE × MUTATION × TREATMENT ───

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            diagnosis_stage="Stage IV",
            treatments_type="Pharmaceutical Therapy, NOS",
            KRAS_protein_change="p.G12D",
            kras_cf=P(lambda x: x is not None),
        ),
    )
    def rule_stage_iv_kras_pharma(self, kras_cf):
        cf_rule = 0.90
        cf_conclusion = cf_rule * kras_cf
        self.declare(
            Result(
                type="Therapy",
                message=" No Response: Stage IV + KRAS G12D + Pharma — linked to aggressive, drug-resistant tumors.",
                certainty_factor=cf_conclusion,
                source="PMID: 35561455",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Consider adding anti-angiogenic therapy or clinical trials for resistant KRAS-driven Stage IV cases.",
                certainty_factor=0.88,
                source="NCCN Guidelines 2024",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            diagnosis_stage="Stage IIA",
            treatments_type="Radiation Therapy, NOS",
            BRAF_status="Mild",
            braf_cf=P(lambda x: x is not None),
        ),
    )
    def rule_stage2_braf(self, braf_cf):
        cf_rule = 0.80
        cf_conclusion = cf_rule * braf_cf
        self.declare(
            Result(
                type="Therapy",
                message=" Response: Early-stage (IIA) with mild BRAF — likely good response to radiotherapy.",
                certainty_factor=cf_conclusion,
                source="PMID: 31739200",
            )
        )

    @Rule(Fact(action="analyze_patient"), Patient(diagnosis_stage="Stage IIA"))
    def rule_stage2_general(self):
        cf_rule = 0.90
        self.declare(
            Result(
                type="Therapy",
                message=" Response: Stage IIA tumors typically have a favorable prognosis.",
                certainty_factor=cf_rule,
                source="NCCN",
            )
        )

    @Rule(Fact(action="analyze_patient"), Patient(diagnosis_stage="Stage IV"))
    def rule_stage4_general(self):
        cf_rule = 0.90
        self.declare(
            Result(
                type="Therapy",
                message=" No Response: Stage IV is often associated with metastasis and poor therapeutic response.",
                certainty_factor=cf_rule,
                source="NCCN",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            diagnosis_stage="Stage IIIC",
            MSI_status=P(lambda x: x >= 1.0),
            msi_cf=MATCH.msi_cf,
        ),
    )
    def rule_stage3_high_msi(self, msi_cf):
        cf_rule = 0.85
        cf_conclusion = cf_rule * msi_cf
        self.declare(
            Result(
                type="Therapy",
                message=" Stable: Stage IIIC with high MSI may benefit from immunotherapy with careful monitoring.",
                certainty_factor=cf_conclusion,
                source="NCCN",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Recommend checkpoint inhibitors and intensive follow-up in Stage IIIC + MSI-high.",
                certainty_factor=0.90,
                source="PMID: 29887279",
            )
        )

    # ─── VITAL STATUS & SURVIVAL ───

    @Rule(
        Fact(action="analyze_patient"),
        Patient(vital_status="Alive", survival_months=P(lambda x: x >= 36)),
    )
    def rule_alive_long_survival(self):
        cf_rule = 0.85
        self.declare(
            Result(
                type="Survival",
                message=" Response: Long survival (≥3 years) indicates treatment success or indolent tumor behavior.",
                certainty_factor=cf_rule,
                source="Oncology Outcome Review",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(vital_status="Dead", survival_months=P(lambda x: x < 6)),
    )
    def rule_short_survival_death(self):
        cf_rule = 0.85
        self.declare(
            Result(
                type="Survival",
                message=" No Response: Death within 6 months suggests rapid progression or therapy failure.",
                certainty_factor=cf_rule,
                source="OncoReview 2023",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(survival_months=P(lambda x: x >= 36), vital_status="Alive"),
    )
    def rule_long_survival(self):
        cf_rule = 0.90
        self.declare(
            Result(
                type="Survival",
                message=" Response: Long survival (≥3 years) is indicative of favorable prognosis or exceptional disease control.",
                certainty_factor=cf_rule,
                source="PMID: 32079319",
            )
        )

    # ─── ADVANCED RULES (IMMUNE + RESPONSE INTERPLAY) ───

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_status="Mild",
            MSI_status=P(lambda x: x < 0.1),
            total_perMB=P(lambda x: x < 2),
            kras_cf=P(lambda x: x is not None),
            msi_cf=MATCH.msi_cf,
        ),
    )
    def rule_low_immune_mild_kras(self, kras_cf, msi_cf):
        cf_rule = 0.75
        cf_conclusion = cf_rule * kras_cf * msi_cf
        self.declare(
            Result(
                type="Immune",
                message=" Stable: Mild KRAS in immunologically cold tumor may offer disease control.",
                certainty_factor=cf_conclusion,
                source="PMID: 35196499",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_status="Mild",
            vital_status="Dead",
            survival_months=P(lambda x: x < 6),
            kras_cf=P(lambda x: x is not None),
        ),
    )
    def rule_conflict_kras_death(self, kras_cf):
        cf_rule = 0.65
        cf_conclusion = cf_rule * kras_cf
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Mild KRAS but early death — review other risk factors or immune evasion. Suggest full clinical/molecular review.",
                certainty_factor=cf_conclusion,
                source="Expert Observation",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            BRAF_status="Severe",
            MSI_status=P(lambda x: x >= 1.5),
            NUMBER_OF_CYCLES=P(lambda x: x >= 12),
            braf_cf=P(lambda x: x is not None),
            msi_cf=MATCH.msi_cf,
        ),
    )
    def rule_braf_high_msi_treatment(self, braf_cf, msi_cf):
        cf_rule = 0.85
        cf_conclusion = cf_rule * braf_cf * msi_cf
        self.declare(
            Result(
                type="Immune",
                message=" Response: Despite severe BRAF, long treatment and high MSI suggest success with immunotherapy.",
                certainty_factor=cf_conclusion,
                source="PMID: 35879849",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_impact="LOW",
            KRAS_status=MATCH.status,
            kras_cf=P(lambda x: x is not None),
        ),
        TEST(lambda status: status in ["Mild", "Moderate"]),
    )
    def rule_kras_low_impact_general(self, kras_cf):
        cf_rule = 0.85
        cf_conclusion = cf_rule * kras_cf
        self.declare(
            Result(
                type="Genetic",
                message=" Response: KRAS mutation with low impact (Mild/Moderate) may retain drug sensitivity.",
                certainty_factor=cf_conclusion,
                source="PMID: 34784796",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            MSI_status=P(lambda x: 0.8 <= x < 1.0), msi_cf=MATCH.msi_cf
        ),
    )
    def rule_moderate_msi(self, msi_cf):
        cf_rule = 0.70
        cf_conclusion = cf_rule * msi_cf
        self.declare(
            Result(
                type="Immune",
                message=" Stable: Moderate MSI level — may provide limited benefit from immune checkpoint inhibitors.",
                certainty_factor=cf_conclusion,
                source="PMID: 34695712",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(NUMBER_OF_CYCLES=P(lambda x: 4 <= x <= 6)),
    )
    def rule_moderate_cycles(self):
        cf_rule = 0.70
        self.declare(
            Result(
                type="Therapy",
                message=" Stable: Moderate treatment exposure (4–6 cycles) may offer disease control — consider extending regimen.",
                certainty_factor=cf_rule,
                source="PMID: 35156114",
            )
        )

    # ─── MODERATE RISK INTERACTIONS ───

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_status="Moderate",
            BRAF_status="Moderate",
            MSI_status=P(lambda x: 0.5 <= x < 1.0),
            total_perMB=P(lambda x: 5.0 <= x < 10.0),
            kras_cf=MATCH.k_cf,
            braf_cf=MATCH.b_cf,
            msi_cf=MATCH.m_cf,
            tmb_cf=MATCH.tmb_cf,
        ),
    )
    def rule_combined_moderate_stability(self, k_cf, b_cf, m_cf, tmb_cf):
        cf_rule = 0.75
        cf_conclusion = cf_rule * k_cf * b_cf * m_cf * tmb_cf
        self.declare(
            Result(
                type="Genetic",
                message=" Stable: Combined moderate mutations and intermediate immune profile suggest partial response potential.",
                certainty_factor=cf_conclusion,
                source="PMID: 35221004",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(BRAF_protein_change="p.K601E"),
        OR(
            Patient(survival_months=P(lambda x: x >= 18)),
            Patient(MSI_status=P(lambda x: x >= 1.0), msi_cf=MATCH.m_cf),
        ),
    )
    def rule_braf_k601e_partial_response(self, m_cf=None):
        cf_rule = 0.80
        cf_conclusion = cf_rule * m_cf if m_cf else cf_rule
        self.declare(
            Result(
                type="Genetic",
                message=" Response: BRAF K601E with favorable survival or MSI may indicate good therapeutic response — consider MEK/BRAF inhibitors.",
                certainty_factor=cf_conclusion,
                source="PMID: 33378690",
            )
        )

    # ─── STABILITY & MODERATE CASES ───


    @Rule(
        Fact(action="analyze_patient"),
        Patient(BRAF_protein_change="p.K601E"),
        Patient(MSI_status=P(lambda x: 0.5 <= x <= 1.0), msi_cf=MATCH.msi_cf),
        Patient(survival_months=P(lambda x: 12 <= x <= 24)),
    )
    def rule_k601e_partial_response(self, msi_cf):
        cf_base = 0.72
        cf_conclusion = cf_base * msi_cf
        self.declare(
            Result(
                type="Stable",
                message=" Stable: K601E variant with moderate MSI and survival suggests partial response — evaluate long-term durability.",
                certainty_factor=cf_conclusion,
                source="PMID: 33378690"
            )
        )


    @Rule(
    Fact(action="analyze_patient"),
    Patient(KRAS_status="Mild"),
    Patient(BRAF_status="Mild"),
    Patient(MSI_status=P(lambda x: x < 0.3), msi_cf=MATCH.msi_cf),
    Patient(total_perMB=P(lambda x: x < 5), tmb_cf=MATCH.tmb_cf),
    Patient(diagnosis_stage=MATCH.stage),
    TEST(lambda stage: stage in ["Stage II", "Stage III"]),
)
    def rule_early_stage_wildtype_low_immune(self, msi_cf, tmb_cf, stage):
        cf_base = 0.68
        cf_conclusion = cf_base * msi_cf * tmb_cf
        self.declare(
            Result(
                type="Stable",
                message=f" Stable: {stage} disease with wild-type profile and low immune burden — recommend standard monitoring every 6–12 months.",
                certainty_factor=cf_conclusion,
                source="PMID: 30865883"
            )
        )
    @Rule(
    Fact(action="analyze_patient"),
    Patient(MSI_status=P(lambda x: 0.4 <= x <= 0.6), msi_cf=MATCH.msi_cf),
    Patient(total_perMB=P(lambda x: 3 <= x <= 10), tmb_cf=MATCH.tmb_cf),
    Patient(survival_months=P(lambda x: x < 24)),
)
    def rule_immune_dormant_window(self, msi_cf, tmb_cf):
        cf_base = 0.65
        cf_conclusion = cf_base * msi_cf * tmb_cf
        self.declare(
            Result(
                type="Stable",
                message=" Stable: Intermediate immune profile with short survival — possible immune dormancy or resistance delay.",
                certainty_factor=cf_conclusion,
                source="PMID: 33140315"
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Recommend immune panel reassessment and repeat MSI profiling.",
                certainty_factor=cf_conclusion + 0.05,
                source="PMID: 33647048"
            )
        )


    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Moderate"),
        Patient(MSI_status=P(lambda x: 0.3 <= x < 0.6)),
        Patient(total_perMB=P(lambda x: 2 <= x < 6)),
        Patient(diagnosis_stage="Stage IIA"),
    )
    def rule_stage2_moderate_stability(self):
        self.declare(
            Result(
                type="Stable",
                message=" Stable: Stage IIA with moderate KRAS and low-intermediate immune markers — disease control likely.",
                certainty_factor=0.72,
                source="NCCN Guidelines 2024",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Moderate"),
        Patient(BRAF_status="Moderate"),
        Patient(MSI_status=P(lambda x: 0.5 <= x < 1.0)),
        Patient(total_perMB=P(lambda x: 5 <= x <= 10)),
        Patient(NUMBER_OF_CYCLES=P(lambda x: 4 <= x <= 6)),
    )
    def rule_stable_combination(self):
        cf_rule = 0.75  # [PMID: 35221004]
        self.declare(
            Result(
                type="Stable",
                message=" Stable: Moderate mutational and immune profile with sufficient therapy suggest disease control.",
                certainty_factor=cf_rule,
                source="PMID: 35221004",
            )
        )

   
    @Rule(
        Fact(action="analyze_patient"),
        Patient(MSI_status=P(lambda x: 0.5 <= x < 1.0)),
    )
    def rule_msi_moderate_response(self):
        cf_rule = 0.70  # [PMID: 34695712]
        self.declare(
            Result(
                type="Immune",
                message=" Stable: Moderate MSI level — may provide limited benefit from immune checkpoint inhibitors.",
                certainty_factor=cf_rule,
                source="PMID: 34695712",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Moderate"),
        Patient(BRAF_status="Moderate"),
        Patient(NUMBER_OF_CYCLES=P(lambda x: 4 <= x <= 6)),
        Patient(MSI_status=P(lambda x: 0.5 <= x <= 1.0)),
        Patient(survival_months=P(lambda x: 12 <= x <= 24)),
    )
    def rule_moderate_all(self):
        cf_rule = 0.75
        self.declare(
            Result(
                type="Therapy",
                message=" Stable: Moderate biomarkers and treatment duration suggest disease stabilization. Continue follow-up every 3–6 months. [PMID: 32507613]",
                certainty_factor=cf_rule,
            )
        )

    # ─── ML-DERIVED RULES (via Decision Tree + SMOTE) ───
    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Moderate"),
        Patient(BRAF_status="Moderate"),
        Patient(MSI_status=P(lambda x: 0.5 <= x <= 0.9)),
        Patient(total_perMB=P(lambda x: 5 <= x <= 10)),
        Patient(NUMBER_OF_CYCLES=P(lambda x: 4 <= x <= 6)),
        Patient(THERAPY_ONGOING="Yes"),
    )
    def rule_ml_partial_response_moderate(self):
        cf_rule = 0.75
        self.declare(
            Result(
                type="ML",
                message=" ML-Predicted: Partial Response — intermediate mutation and immune profile suggest partial control under ongoing therapy.",
                certainty_factor=cf_rule,
                source="ML Pattern | NCCN-Calibrated Thresholds",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(THERAPY_ONGOING="No"),
        Patient(BRAF_impact="LOW"),
        Patient(MSI_status=P(lambda x: x <= 0.62)),
        Patient(diagnosis_stage="Stage IIA"),
    )
    def rule_ml_complete_response(self):
        cf_rule = 0.85  # From ML model (high confidence)
        self.declare(
            Result(
                type="ML",
                message=" ML-Predicted: Complete Response — Early stage, low MSI & BRAF impact align with favorable outcome. [Source: Tree-based model]",
                certainty_factor=cf_rule,
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(THERAPY_ONGOING="No"),
        Patient(BRAF_impact="LOW"),
        Patient(MSI_status=P(lambda x: x <= 0.62)),
        Patient(diagnosis_stage="Stage IIIC"),
    )
    def rule_ml_partial_response(self):
        cf_rule = 0.75  # ML model prediction
        self.declare(
            Result(
                type="ML",
                message=" ML-Predicted: Partial Response — Advanced stage despite favorable BRAF/MSI may reduce efficacy.",
                certainty_factor=cf_rule,
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(THERAPY_ONGOING="No"),
        Patient(BRAF_impact="LOW"),
        Patient(MSI_status=P(lambda x: x > 0.62)),
        Patient(survival_months=P(lambda x: x <= 11.99)),
    )
    def rule_ml_clinical_progressive(self):
        cf_rule = 0.80  # ML model – progression with biomarker mismatch
        self.declare(
            Result(
                type="ML",
                message=" ML-Predicted: Clinical Progressive Disease — High MSI with short survival suggests progression.",
                certainty_factor=cf_rule,
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(THERAPY_ONGOING="Yes"),
        Patient(KRAS_status="Mild"),
        Patient(total_perMB=P(lambda x: x > 6)),
        Patient(patient_id=P(lambda x: x <= 292.5)),
    )
    def rule_ml_radiographic_progressive(self):
        cf_rule = 0.80  # ML rule
        self.declare(
            Result(
                type="ML",
                message=" ML-Predicted: Radiographic Progressive Disease — High TMB + therapy may indicate subclinical progression.",
                certainty_factor=cf_rule,
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Moderate"),
        Patient(BRAF_status="Mild"),
        Patient(survival_months=P(lambda x: x >= 36)),
    )
    def rule_stable_response_moderate_kras(self):
        cf_rule = 0.75
        self.declare(
            Result(
                type="Survival",
                message=" Stable: Moderate KRAS with mild BRAF and long survival suggests durable disease control.",
                certainty_factor=cf_rule,
            )
        )

    # ─── ADVANCED ANALYSIS – ENHANCEMENTS ───

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_protein_change="p.G12D",
            KRAS_status="Mild",
            KRAS_impact="LOW",
            vital_status="Alive",
            survival_months=P(lambda x: x >= 36),
            kras_cf=MATCH.k_cf,
        ),
    )
    def rule_conflict_kras_g12d_mild_but_good(self, k_cf):
        cf_rule = 0.65
        cf_conclusion = cf_rule * k_cf
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: KRAS G12D is usually resistant, but mild status and long survival suggest favorable subclone or treatment effect.",
                certainty_factor=cf_conclusion,
                source="PMID: 36846823",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Consider tumor heterogeneity or epigenetic modulation.",
                certainty_factor=0.80,
                source="PMID: 36846823",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_status="Severe",
            KRAS_impact="HIGH",
            vital_status="Alive",
            survival_months=P(lambda x: x >= 36),
            kras_cf=MATCH.k_cf,
        ),
    )
    def rule_conflict_kras_severe_but_survived(self, k_cf):
        cf_rule = 0.65
        cf_conclusion = cf_rule * k_cf
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Severe KRAS mutation but long survival — possible exception due to unknown modifiers or immune activity.",
                certainty_factor=cf_conclusion,
                source="PMID: 32079319",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggest immune checkpoint evaluation or retrospective treatment mapping.",
                certainty_factor=0.78,
                source="PMID: 32079319",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_status="Moderate",
            BRAF_status="Moderate",
            MSI_status=P(lambda x: 0.5 <= x <= 1.0),
            total_perMB=P(lambda x: 3 <= x <= 10),
            kras_cf=MATCH.k_cf,
            braf_cf=MATCH.b_cf,
            msi_cf=MATCH.m_cf,
            tmb_cf=MATCH.t_cf,
        ),
    )

    # ─── EXCEPTIONAL SURVIVAL CASES ───
    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Mild"),
        Patient(BRAF_status="Mild"),
        Patient(MSI_status=P(lambda x: x < 0.3)),
        Patient(total_perMB=P(lambda x: x < 2)),
        Patient(survival_months=P(lambda x: x >= 80)),
    )
    def rule_exceptional_low_mutation_but_long_survival(self):
        self.declare(
            Result(
                type="Exceptional",
                message=" Exceptional: Mild mutations, low MSI/TMB, yet long survival — consider epigenetic regulation or slow-growing tumor subtype.",
                certainty_factor=0.75,
                source="PMID: 34378912",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(diagnosis_stage="Stage IV"),
        Patient(survival_months=P(lambda x: x >= 36)),
    )
    def rule_stage4_long_survival(self):
        cf_rule = 0.70
        self.declare(
            Result(
                type="Exceptional",
                message=" Exceptional Case: Stage IV disease with long survival (≥36 months) — consider molecular review or immune profiling. [PMID: 35788316]",
                certainty_factor=cf_rule,
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(diagnosis_stage="Stage IIIC"),
        Patient(survival_months=P(lambda x: x >= 36)),
    )
    def rule_stage_iiic_long_survival(self):
        cf_rule = 0.70
        self.declare(
            Result(
                type="Exceptional",
                message=" Rare Case: Stage IIIC with survival beyond expectations — explore immune microenvironment and mutation burden. [PMID: 32492650]",
                certainty_factor=cf_rule,
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Severe"),
        Patient(KRAS_impact="HIGH"),
        Patient(BRAF_status="Severe"),
        Patient(BRAF_impact="HIGH"),
        Patient(MSI_status=P(lambda x: x < 0.3)),
        Patient(total_perMB=P(lambda x: x < 5)),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 4)),
        Patient(survival_months=P(lambda x: x >= 120)),
        Patient(vital_status="Alive"),
    )
    def rule_extreme_conflict_kras_braf_survival(self):
        cf_rule = 0.60
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Dual high-risk mutations (KRAS & BRAF), low immune markers, short treatment — yet long survival.",
                certainty_factor=cf_rule,
            )
        )
        self.declare(
            Result(
                type="Exceptional",
                message=" Note: Possible immunologic or rare genetic modifier protecting patient — recommend deep molecular profiling. [PMID: 34710943]",
                certainty_factor=cf_rule,
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggest: Archive for exceptional response research, consider long-term immune memory role.",
                certainty_factor=0.90,
                source="PMID: 37059545",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Severe"),
        Patient(BRAF_status="Severe"),
        Patient(MSI_status=P(lambda x: x < 0.5)),
        Patient(total_perMB=P(lambda x: x < 5)),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 4)),
        Patient(survival_months=P(lambda x: x >= 120)),
        Patient(THERAPY_ONGOING="No"),
    )
    def rule_extreme_exceptional_survival(self):
        cf_rule = 0.60
        self.declare(
            Result(
                type="Exceptional",
                message=" Exceptional: Severe dual mutation, low immune markers, short therapy — yet survival >10y suggests rare resistance or genomic uniqueness.",
                certainty_factor=cf_rule,
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggest archiving and full profiling: germline variants, immunogenetics, and retrospective clinical audit.",
                certainty_factor=0.88,
                source="PMID: 37059545",
            )
        )

    # ─── CONFLICTS & EXCEPTIONAL RESPONSES ───
    @Rule(
    Fact(action="analyze_patient"),
    Patient(KRAS_status="Moderate"),
    Patient(MSI_status=P(lambda x: x < 0.5), msi_cf=MATCH.msi_cf),
    Patient(total_perMB=P(lambda x: x < 5), tmb_cf=MATCH.tmb_cf),
    Patient(survival_months=P(lambda x: x >= 60)),
)
    def rule_kras_moderate_but_long_survival(self, msi_cf, tmb_cf):
        cf_base = 0.65
        cf_conclusion = cf_base * msi_cf * tmb_cf
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: KRAS moderate with low immune markers but long survival — consider late immune rebound or unmeasured modifiers.",
                certainty_factor=cf_conclusion,
                source="PMID: 35910014"
            )
        )
    
    @Rule(
    Fact(action="analyze_patient"),
    Patient(diagnosis_stage="Stage IIIC"),
    Patient(KRAS_status="Mild"),
    Patient(BRAF_status="Mild"),
    Patient(MSI_status=P(lambda x: x < 0.4), msi_cf=MATCH.msi_cf),
    Patient(total_perMB=P(lambda x: x < 5), tmb_cf=MATCH.tmb_cf),
    Patient(survival_months=P(lambda x: x >= 48)),
)
    def rule_iiic_long_survival_wild(self, msi_cf, tmb_cf):
        cf_base = 0.68
        cf_conclusion = cf_base * msi_cf * tmb_cf
        self.declare(
            Result(
                type="Exceptional",
                message=" Exceptional: Stage IIIC with minimal mutation burden and long survival — spontaneous control or underestimated immune response.",
                certainty_factor=cf_conclusion,
                source="PMID: 32492650"
            )
        )


    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_protein_change="p.G12D",
            KRAS_status="Mild",
            survival_months=P(lambda x: x >= 36),
            kras_cf=MATCH.k_cf,
        ),
    )
    def rule_conflict_kras_g12d_mild_survival(self, k_cf):
        cf_rule = 0.65
        cf_conclusion = cf_rule * k_cf
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: KRAS G12D is usually resistant, but mild status and long survival suggest an exception.",
                certainty_factor=cf_conclusion,
                source="PMID: 36846823",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Consider: tumor heterogeneity or subclonal KRAS behavior.",
                certainty_factor=0.80,
                source="PMID: 36846823",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_status="Severe",
            KRAS_impact="HIGH",
            survival_months=P(lambda x: x >= 36),
            kras_cf=MATCH.k_cf,
        ),
    )
    def rule_conflict_severe_kras_but_survives(self, k_cf):
        cf_rule = 0.65
        cf_conclusion = cf_rule * k_cf
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Severe KRAS mutation but long survival — potential treatment success despite poor markers.",
                certainty_factor=cf_conclusion,
                source="PMID: 32079319",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Recommend molecular re-evaluation or immunotherapy exposure history.",
                certainty_factor=0.75,
                source="PMID: 32079319",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_status="Severe",
            KRAS_protein_change="p.G12D",
            NUMBER_OF_CYCLES=P(lambda x: x <= 3),
            survival_months=P(lambda x: x >= 48),
            kras_cf=MATCH.k_cf,
        ),
    )
    def rule_paradoxical_survival_with_risk(self, k_cf):
        cf_rule = 0.60
        cf_conclusion = cf_rule * k_cf
        self.declare(
            Result(
                type="Exceptional",
                message=" Note: Despite high-risk markers (Severe KRAS, G12D, short treatment), long survival suggests rare resistance reversal.",
                certainty_factor=cf_conclusion,
                source="PMID: 37059545",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggested Action: Archive this case for review or potential case study.",
                certainty_factor=0.85,
                source="PMID: 37059545",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            survival_months=P(lambda x: x >= 48),
            KRAS_status="Severe",
            THERAPY_ONGOING="No",
            kras_cf=MATCH.k_cf,
        ),
    )
    def rule_exceptional_survival(self, k_cf):
        cf_rule = 0.65
        cf_conclusion = cf_rule * k_cf
        self.declare(
            Result(
                type="Exceptional",
                message=" Exceptional Case: Severe mutation with early therapy stop but long survival — consider immune or genetic resilience.",
                certainty_factor=cf_conclusion,
                source="PMID: 35788316",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggest: Genomic sequencing or immune profiling.",
                certainty_factor=0.88,
                source="PMID: 35788316",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_status="Severe",
            survival_months=P(lambda x: x >= 48),
            kras_cf=MATCH.k_cf,
        ),
    )
    def rule_merged_conflict_kras_survival(self, k_cf):
        cf_rule = 0.60
        cf_conclusion = cf_rule * k_cf
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Severe KRAS mutation with unexpected long survival — resistance may be overridden.",
                certainty_factor=cf_conclusion,
                source="PMID: 36950031",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Review for immune escape or off-target benefit from combination therapy.",
                certainty_factor=0.80,
                source="PMID: 36950031",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            KRAS_status="Severe",
            survival_months=P(lambda x: x >= 36),
            kras_cf=MATCH.k_cf,
        ),
    )
    def rule_conflict_kras_survival(self, k_cf):
        cf_rule = 0.65
        cf_conclusion = cf_rule * k_cf
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Severe KRAS mutation but long survival — consider resistance reversal or host immune factors.",
                certainty_factor=cf_conclusion,
                source="PMID: 36950031",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(MSI_status=P(lambda x: x < 0.5)),
        Patient(msi_cf=P(lambda x: x > 0.9)),
        Patient(survival_months=P(lambda x: x > 36)),
    )
    def rule_conflict_msi_confidence_vs_outcome(self):
        cf_rule = 0.60
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: MSI reported as low with high input confidence, but survival is long — review MSI methodology or patient-specific immunity.",
                certainty_factor=cf_rule,
                source="PMID: 33140315",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Severe"),
        Patient(BRAF_status="Severe"),
        Patient(MSI_status=P(lambda x: x < 0.5)),
        Patient(total_perMB=P(lambda x: x < 5)),
    )
    def rule_conflict_low_immune_and_severe(self):
        cf_rule = 0.70
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Severe KRAS/BRAF with low MSI/TMB — likely immune evasion and high resistance.",
                certainty_factor=cf_rule,
                source="PMID: 29887279",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(MSI_status=P(lambda x: x < 0.3)),
        Patient(total_perMB=P(lambda x: x < 5)),
        Patient(survival_months=P(lambda x: x >= 100)),
    )
    def rule_exceptional_low_immune_high_survival(self):
        cf_rule = 0.75
        self.declare(
            Result(
                type="Exceptional",
                message=" Exceptional: Long survival despite low immune markers — suggests unique resistance or immune-independent response.",
                certainty_factor=cf_rule,
                source="PMID: 35788316",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Severe"),
        Patient(MSI_status=P(lambda x: x < 0.3)),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 3)),
    )
    def rule_high_risk_cluster(self):
        cf_rule = 0.80
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Severe KRAS + low immune score + short therapy — high risk of recurrence or progression.",
                certainty_factor=cf_rule,
                source="PMID: 36950031",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(THERAPY_ONGOING="No"),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 4)),
        Patient(survival_months=P(lambda x: x >= 36)),
    )
    def rule_paradoxical_survival(self):
        cf_rule = 0.60
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Early therapy stop but long survival — reassess potential immune-mediated control or patient-specific factors.",
                certainty_factor=cf_rule,
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggest immune profiling or retrospective review for exceptional outcome.",
                certainty_factor=0.80,
                source="PMID: 34237709",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(survival_months=P(lambda x: x < 6)),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x >= 6)),
    )
    def rule_early_death_despite_treatment(self):
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Rapid disease progression despite extended treatment — possible intrinsic resistance.",
                certainty_factor=0.70,  # [PMID: 32345678]
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggest case review for hyperprogression or early escape.",
                certainty_factor=0.80,
                source="PMID: 32345678",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(BRAF_protein_change="p.V600M"),
        Patient(survival_months=P(lambda x: x < 6)),
    )
    def rule_rare_braf_v600m_poor_outcome(self):
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Rare BRAF V600M with early death — prognosis may be worse than V600E variants.",
                certainty_factor=0.65,
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Recommend: Investigate non-V600E BRAF variants for altered therapy resistance.",
                certainty_factor=0.75,
                source="PMID: 28167255",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(diagnosis_stage="Stage IV"),
        Patient(treatments_type=MATCH.treat),
        TEST(lambda treat: "Radiation" in treat),
        Patient(survival_months=P(lambda x: x >= 18)),
    )
    def rule_stage4_radiation_survival(self):
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Prolonged survival in Stage IV with radiation only — review for immune resilience or minimal disease burden.",
                certainty_factor=0.70,
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggest PET-CT reassessment and consider immune/genomic contributors.",
                certainty_factor=0.75,
                source="PMID: 31928476",
            )
        )

    # ─── ADVANCED COMPOSITE RULES ───

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Severe", KRAS_impact="HIGH", kras_cf=MATCH.k_cf),
        Patient(MSI_status=P(lambda x: x >= 1.0), msi_cf=MATCH.m_cf),
        Patient(total_perMB=P(lambda x: x >= 20), tmb_cf=MATCH.t_cf),
    )
    def rule_kras_severe_but_high_immune(self, k_cf, m_cf, t_cf):
        cf_rule = 0.75
        cf_conclusion = cf_rule * k_cf * m_cf * t_cf
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Despite severe KRAS, high MSI and TMB may override resistance.",
                certainty_factor=cf_conclusion,
                source="PMID: 29887279",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Consider prioritizing immunotherapy or combined approaches.",
                certainty_factor=0.85,
                source="FDA, NCCN",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        OR(
            Patient(MSI_status=P(lambda x: x >= 10), msi_cf=MATCH.m_cf),
            Patient(total_perMB=P(lambda x: x >= 30), tmb_cf=MATCH.t_cf),
        ),
        Patient(THERAPY_ONGOING="No"),
    )
    def rule_high_immune_either(self, m_cf=None, tmb_cf=None):
        cf_rule = 0.85
        cf_conclusion = cf_rule * (m_cf or 1.0) * (tmb_cf or 1.0)
        self.declare(
            Result(
                type="Recommendation",
                message=" Strong immune biomarker (MSI or TMB) detected — initiate immunotherapy or re-challenge with anti-PD-1.",
                certainty_factor=cf_conclusion,
                source="PMID: 29887279",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status=MATCH.k),
        Patient(BRAF_status=MATCH.b),
        Patient(survival_months=P(lambda x: x >= 24)),
        TEST(lambda k, b: k in ["Moderate", "Severe"] and b in ["Moderate", "Severe"]),
    )
    def rule_stabilized_disease_despite_mutations(self):
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Moderate/Severe mutations but long survival — indicates partial or off-target control.",
                certainty_factor=0.65,
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Consider retrospective mapping of therapy combinations or resistance modulation.",
                certainty_factor=0.80,
                source="PMID: 37059545",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(THERAPY_ONGOING="No"),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 4)),
        OR(
            Patient(MSI_status=P(lambda x: x >= 15)),
            Patient(total_perMB=P(lambda x: x >= 30)),
        ),
        Patient(survival_months=P(lambda x: x >= 24)),
    )
    def rule_early_stop_but_success_immune(self):
        self.declare(
            Result(
                type="Exceptional",
                message=" Note: Early stop but high immune biomarkers and long survival — likely immune-mediated control.",
                certainty_factor=0.80,
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Recommend follow-up with immune profiling and tumor escape mechanism screening.",
                certainty_factor=0.85,
                source="PMID: 33140315",
            )
        )

    # ─── IMPROVED CONFLICT DETECTION ───
    @Rule(
    Fact(action="analyze_patient"),
    OR(
        Patient(MSI_status=P(lambda x: x >= 10), msi_cf=MATCH.mcf),
        Patient(total_perMB=P(lambda x: x >= 30), tmb_cf=MATCH.tcf),
    ),
    Patient(THERAPY_ONGOING="No"),
    Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 4)),
    Patient(survival_months=P(lambda x: x >= 36)),
)
    def rule_immunogenic_control(self, mcf=1.0, tcf=1.0):
        base_cf = 0.78
        cf_conclusion = base_cf * max(mcf, tcf)
        self.declare(
            Result(
                type="Exceptional",
                message=" Exceptional: Very high immune signature with short therapy but long survival — likely immune-driven tumor control.",
                certainty_factor=cf_conclusion,
                source="PMID: 33140315"
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Recommend: Archive for immune memory profiling and long-term monitoring.",
                certainty_factor=cf_conclusion + 0.07,
                source="FDA, 2020 Approval"
            )
        )
    @Rule(
    Fact(action="analyze_patient"),
    Patient(KRAS_status="Severe"),
    OR(
        Patient(BRAF_status="Moderate"),
        Patient(MSI_status=P(lambda x: x > 0.9), msi_cf=MATCH.mcf),
    ),
    Patient(survival_months=P(lambda x: x >= 36)),
)
    def rule_kras_severe_partial_modifiers(self, mcf=1.0):
        base_cf = 0.66
        cf_conclusion = base_cf * mcf
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Severe KRAS with partial modifying factors (BRAF moderate or MSI high) and extended survival.",
                certainty_factor=cf_conclusion,
                source="PMID: 35910014"
            )
        )
    @Rule(
        Fact(action="analyze_patient"),
        Patient(diagnosis_stage=MATCH.stage),
        TEST(lambda stage: stage in ["Stage I", "Stage II"]),
        OR(
            Patient(BRAF_protein_change="p.V600E"),
            Patient(BRAF_protein_change="p.K601E"),
        ),
        AND(
            Patient(MSI_status=P(lambda x: x >= 0.6), msi_cf=MATCH.mcf),
            Patient(total_perMB=P(lambda x: x >= 10), tmb_cf=MATCH.tcf),
        ),
    )
    def rule_early_stage_braf_and_immune(self, mcf, tcf):
        base_cf = 0.77
        cf_conclusion = base_cf * mcf * tcf
        self.declare(
            Result(
                type="Stable",
                message=" Stable: Early stage BRAF-mutant with moderate-high immune markers — potential for durable response.",
                certainty_factor=cf_conclusion,
                source="NCCN 2024"
            )
        )
    @Rule(
    Fact(action="analyze_patient"),
    Patient(MSI_status=P(lambda x: x < 0.3), msi_cf=MATCH.mcf),
    OR(
        Patient(THERAPY_ONGOING="Yes"),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 4)),
    ),
    Patient(survival_months=P(lambda x: x >= 60)),
)
    def rule_silent_immunity_or_host_factors(self, mcf):
        base_cf = 0.70
        cf_conclusion = base_cf * mcf
        self.declare(
            Result(
                type="Exceptional",
                message=" Exceptional: Low MSI with ongoing therapy or short therapy but long survival — may reflect host immunity or tumor dormancy.",
                certainty_factor=cf_conclusion,
                source="PMID: 37059545"
            )
        )
    @Rule(
    Fact(action="analyze_patient"),
    OR(
        Patient(BRAF_status="Moderate"),
        Patient(MSI_status=P(lambda x: x >= 0.7), msi_cf=MATCH.mcf),
    ),
    Patient(survival_months=P(lambda x: 24 <= x <= 36)),
)
    def rule_partial_response_braf_or_msi(self, mcf=1.0):
        base_cf = 0.65
        cf_conclusion = base_cf * mcf
        self.declare(
            Result(
                type="Stable",
                message=" Stable: Intermediate survival with partial immune/genetic signals (MSI or BRAF) — consider durable disease control.",
                certainty_factor=cf_conclusion,
                source="PMID: 33378690"
            )
        )
        
    
    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Severe"),
        Patient(KRAS_impact="HIGH"),
        Patient(survival_months=P(lambda x: x >= 36)),
    )
    def rule_conflict_severe_kras_long_survival(self):
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Severe KRAS mutation but long survival — may suggest resistance override or exceptional response.",
                certainty_factor=0.70,
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Recommend further analysis of tumor heterogeneity or immune compensation.",
                certainty_factor=0.75,
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_protein_change="p.G12D"),
        Patient(KRAS_status=MATCH.status),
        TEST(lambda status: status in ["Mild", "Moderate"]),
        Patient(survival_months=P(lambda x: x >= 36)),
    )
    def rule_conflict_kras_g12d_but_long_survival(self):
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: KRAS G12D mutation is typically resistant, but long survival suggests potential subclonal behavior.",
                certainty_factor=0.70,
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(THERAPY_ONGOING="No"),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 4)),
        Patient(survival_months=P(lambda x: x >= 48)),
    )
    def rule_conflict_short_therapy_long_survival(self):
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Short treatment course but extended survival — reassess potential immune or genomic factors.",
                certainty_factor=0.70,
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggest immune profiling or germline sequencing.",
                certainty_factor=0.80,
                source="PMID: 35788316",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Severe"),
        Patient(THERAPY_ONGOING="No"),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 4)),
        Patient(survival_months=P(lambda x: x >= 48)),
    )
    def rule_conflict_combined_kras_therapy_survival(self):
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Severe KRAS + short therapy + long survival — highly atypical outcome.",
                certainty_factor=0.65,
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggest archiving this case for in-depth review or case report.",
                certainty_factor=0.80,
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(MSI_status=P(lambda x: x < 0.3), msi_cf=MATCH.m_cf),
        Patient(total_perMB=P(lambda x: x < 2), tmb_cf=MATCH.t_cf),
        Patient(survival_months=P(lambda x: x >= 36)),
    )
    def rule_conflict_low_immune_markers_but_survives(self, m_cf, t_cf):
        cf_rule = 0.60
        cf_conclusion = cf_rule * m_cf * t_cf
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Low MSI/TMB but prolonged survival — contradicts expected immune responsiveness.",
                certainty_factor=cf_conclusion,
                source="PMID: 29887279",
            )
        )
        self.declare(
            Result(
                type="Recommendation",
                message="🔎 Consider immune evasion resistance or misclassification.",
                certainty_factor=0.70,
                source="PMID: 29887279",
            )
        )

    # ─── THERAPEUTIC RECOMMENDATIONS ───
    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Moderate"),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 4)),
        Patient(diagnosis_stage="Stage IIA"),
        Patient(THERAPY_ONGOING="Yes"),
    )
    def rule_suggest_extend_chemo_stage2(self):
        self.declare(
            Result(
                type="Recommendation",
                message=" Stage IIA with moderate KRAS and low cycles — consider extending chemotherapy course to ≥6 for optimal benefit.",
                certainty_factor=0.78,
                source="PMID: 35156114",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Moderate"),
        Patient(BRAF_status="Moderate"),
        Patient(MSI_status=P(lambda x: 0.5 <= x <= 1.0)),
        Patient(survival_months=P(lambda x: 12 <= x <= 24)),
    )
    def rule_recommend_follow_up_moderate_case(self):
        cf_rule = 0.70  # [NCCN 2024]
        self.declare(
            Result(
                type="Recommendation",
                message=" Recommend: Continue 3–6 month monitoring for moderate-risk profiles with stable course.",
                certainty_factor=cf_rule,
                source="NCCN Guidelines 2024",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(BRAF_protein_change="p.K601E"),
        Patient(diagnosis_stage="Stage IIIC"),
        Patient(MSI_status=P(lambda x: 0.5 <= x <= 1.0)),
        Patient(THERAPY_ONGOING="Yes"),
    )
    def rule_k601e_radiation_combo(self):
        cf_rule = 0.75  # [PMID: 35879849]
        self.declare(
            Result(
                type="Recommendation",
                message=" BRAF K601E in Stage IIIC + radiation may benefit from MEK inhibitor addition — consider multimodal strategy.",
                certainty_factor=cf_rule,
                source="PMID: 35879849",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(
            MSI_status=P(lambda x: x < 0.5),
            total_perMB=P(lambda x: x < 5),
            msi_cf=MATCH.m_cf,
            tmb_cf=MATCH.t_cf,
        ),
        Patient(treatments_type="Pharmaceutical Therapy, NOS"),
    )
    def suggest_add_immunotherapy(self, m_cf, t_cf):
        cf_rule = 0.75
        cf_conclusion = cf_rule * m_cf * t_cf
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggestion: Low MSI/TMB with pharma-only therapy — consider immunotherapy addition (e.g. checkpoint inhibitors).",
                certainty_factor=cf_conclusion,
                source="PMID: 31483745",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(BRAF_protein_change="p.V600E"),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x < 6)),
    )
    def suggest_extend_braf(self):
        cf_rule = 0.80
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggestion: V600E mutation detected — consider extending BRAF-targeted treatment.",
                certainty_factor=cf_rule,
                source="PMID: 35221741",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(diagnosis_stage="Stage IIIC"),
        Patient(MSI_status=P(lambda x: x >= 1.0), msi_cf=MATCH.m_cf),
    )
    def suggest_monitor_stage_iiic(self, m_cf):
        cf_rule = 0.85
        cf_conclusion = cf_rule * m_cf
        self.declare(
            Result(
                type="Recommendation",
                message=" Recommendation: Stage IIIC + High MSI — initiate intensified monitoring or combine immunotherapy.",
                certainty_factor=cf_conclusion,
                source="NCCN Guidelines 2024",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x <= 3)),
        Patient(survival_months=P(lambda x: x >= 36)),
    )
    def rule_short_therapy_long_survival(self):
        self.declare(
            Result(
                type="Exceptional",
                message=" Note: Long survival despite very short therapy — may reflect immune-mediated control or tumor indolence.",
                certainty_factor=0.80,
                source="PMID: 34237709",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x >= 10)),
        Patient(survival_months=P(lambda x: x < 24)),
    )
    def rule_long_therapy_low_survival(self):
        self.declare(
            Result(
                type="Conflict",
                message=" Conflict: Extended therapy without proportional survival benefit — may indicate drug resistance or treatment inefficacy.",
                certainty_factor=0.65,
                source="PMID: 32079319",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(treatments_type=MATCH.t),
        TEST(lambda t: "Chemotherapy" in t),
        Patient(BRAF_protein_change="p.V600E"),
    )
    def rule_combine_braf_chemo(self):
        self.declare(
            Result(
                type="Recommendation",
                message=" BRAF V600E detected — consider combining chemotherapy with BRAF-targeted agents.",
                certainty_factor=0.85,
                source="PMID: 29996028, NCCN",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(MSI_status=P(lambda x: x >= 1.0), msi_cf=MATCH.cf_user),
    )
    def suggest_msi_high_immunotherapy(self, cf_user):
        cf_rule = 0.90
        self.declare(
            Result(
                type="Recommendation",
                message=" Consider immunotherapy (e.g., anti-PD-1 like pembrolizumab) in MSI-high patients.",
                certainty_factor=cf_rule * cf_user,
                source="PMID: 29887279, NCCN Guidelines 2024",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(total_perMB=P(lambda x: x > 20), tmb_cf=MATCH.cf_user),
    )
    def suggest_tmb_high_immunotherapy(self, cf_user):
        cf_rule = 0.90
        self.declare(
            Result(
                type="Recommendation",
                message=" Consider immune checkpoint blockade in high TMB tumors (>20 mut/Mb).",
                certainty_factor=cf_rule * cf_user,
                source="FDA Approval: pembrolizumab, 2020",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(BRAF_protein_change="p.V600E"),
    )
    def suggest_braf_v600e_targeting(self):
        self.declare(
            Result(
                type="Recommendation",
                message=" Use BRAF inhibitors (e.g., vemurafenib, encorafenib) for BRAF V600E mutation.",
                certainty_factor=0.88,
                source="PMID: 26559394, NCCN",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(BRAF_protein_change="p.V600E"),
        Patient(NUMBER_OF_CYCLES=P(lambda x: x < 6)),
    )
    def suggest_extend_braf_therapy(self):
        self.declare(
            Result(
                type="Recommendation",
                message=" Suggest extending duration of BRAF-targeted therapy for incomplete exposure.",
                certainty_factor=0.80,
                source="PMID: 35221741",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(BRAF_protein_change="p.K601E"),
        OR(
            Patient(survival_months=P(lambda x: x >= 18)),
            Patient(MSI_status=P(lambda x: x >= 1.0)),
        ),
    )
    def suggest_braf_k601e_option(self):
        self.declare(
            Result(
                type="Recommendation",
                message=" BRAF K601E with favorable outcome — consider MEK inhibitor-based strategies.",
                certainty_factor=0.75,
                source="PMID: 33378690",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(MSI_status=P(lambda x: x < 0.5), msi_cf=MATCH.msi_cf),
        Patient(total_perMB=P(lambda x: x < 5), tmb_cf=MATCH.tmb_cf),
        Patient(treatments_type="Pharmaceutical Therapy, NOS"),
    )
    def suggest_add_immunotherapy_low_msi_tmb(self, msi_cf, tmb_cf):
        cf_rule = 0.70
        self.declare(
            Result(
                type="Recommendation",
                message=" Low MSI/TMB with pharma therapy — consider adding checkpoint inhibitors.",
                certainty_factor=cf_rule * msi_cf * tmb_cf,
                source="PMID: 31483745",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(diagnosis_stage="Stage IIIC"),
        Patient(MSI_status=P(lambda x: x >= 1.0)),
    )
    def suggest_advanced_msi_monitoring(self):
        self.declare(
            Result(
                type="Recommendation",
                message=" Stage IIIC + MSI-high: Initiate aggressive monitoring + immunotherapy consideration.",
                certainty_factor=0.85,
                source="NCCN Guidelines",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Severe"),
        Patient(KRAS_protein_change="p.G12D"),
    )
    def suggest_anti_egfr_ineffective(self):
        self.declare(
            Result(
                type="Recommendation",
                message=" Do not use anti-EGFR therapy in KRAS G12D or severe KRAS mutations.",
                certainty_factor=0.90,
                source="OncoKB, NCCN",
            )
        )

    # ─── TARGETED THERAPY & DURATION RULES ───

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_protein_change=MATCH.kp),
        TEST(lambda kp: kp in ["p.G12D", "p.G12V"]),
    )
    def rule_anti_egfr_contraindicated(self):
        self.declare(
            Result(
                type="Recommendation",
                message=" Avoid anti-EGFR therapy in KRAS G12D or G12V — known resistance.",
                certainty_factor=0.90,
                source="PMID: 35521543, OncoKB",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(KRAS_status="Wild-Type"),
    )
    def rule_egfr_suitable(self):
        self.declare(
            Result(
                type="Recommendation",
                message=" Consider EGFR inhibitors (e.g., cetuximab) in KRAS wild-type patients.",
                certainty_factor=0.85,
                source="NCCN 2024",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(BRAF_protein_change="p.V600E"),
    )
    def rule_braf_mek_combination(self):
        self.declare(
            Result(
                type="Recommendation",
                message=" BRAF V600E detected — combine BRAF and MEK inhibitors (e.g., encorafenib + binimetinib).",
                certainty_factor=0.88,
                source="PMID: 26559394, NCCN",
            )
        )

    @Rule(
        Fact(action="analyze_patient"),
        Patient(BRAF_protein_change="p.K601E"),
    )
    def rule_braf_k601e_consideration(self):
        self.declare(
            Result(
                type="Recommendation",
                message=" BRAF K601E mutation — consider MEK inhibitor sensitivity evaluation.",
                certainty_factor=0.75,
                source="PMID: 33378690",
            )
        )

        # Exceptional Cases
    @Rule(
    Fact(action="analyze_patient"),
    Patient(MSI_status=P(lambda x: x < 0.3)),
    Patient(total_perMB=P(lambda x: x < 5)),
    Patient(survival_months=P(lambda x: x >= 60)),
    Patient(THERAPY_ONGOING="No"),
)
    def rule_immune_independent_mechanism(self):
        cf = 0.80
        self.declare(Result(
            type="Exceptional",
            message=" Exceptional: Long-term survival with low immune profile and no therapy — possible immune dormancy or alternate resistance.",
            certainty_factor=cf,
            source="PMID: 36515490"
        ))
    @Rule(
    Fact(action="analyze_patient"),
    Patient(diagnosis_stage=MATCH.stage),
    TEST(lambda stage: stage in ["Stage II", "Stage III", "Stage IIA", "Stage IIIC"]),
    Patient(KRAS_status="Mild"),
    Patient(BRAF_status="Mild"),
    Patient(MSI_status=P(lambda x: x >= 0.5)),
    Patient(NUMBER_OF_CYCLES=P(lambda x: x >= 6)),
)
    def rule_early_stage_wildtype_response(self):
        cf = 0.85
        self.declare(Result(
            type="Stable",
            message=" Stage II/III wild-type tumors with immune signal and full therapy often show durable response.",
            certainty_factor=cf,
            source="NCCN Guidelines 2024"
        ))

    @Rule(
        Fact(action="analyze_patient"),
        Patient(MSI_status=P(lambda x: x < 0.5)),
        Patient(total_perMB=P(lambda x: x < 5)),
        Patient(survival_months=P(lambda x: x > 100)),
    )
    def rule_extreme_survival_low_immune(self):
        cf_rule = 0.75
        self.declare(
            Result(
                type="Exceptional",
                message=" Exceptional: Long survival despite low immune markers — suggests unique resistance or immune-independent response.",
                certainty_factor=cf_rule,
                source="PMID: 36950031",
            )
        )

    # ─── DEFAULT ───

    @Rule(Fact(action="analyze_patient"))
    def rule_fallback(self):
        if not any(self.outputs[key] for key in self.outputs if key != "Fallback"):
            self.outputs["Fallback"].append(
                " No specific rule matched — recommend full clinical review with multidisciplinary team."
            )

    def print_final_report(self):
        patient_fact = next((f for f in self.facts.values() if isinstance(f, Patient)), {})
        patient_id = patient_fact.get("patient_id", None)

        try:
            patient_number = int(float(patient_id)) if patient_id else "Unknown"
        except ValueError:
            patient_number = "Unknown"

        print(f"\n{'='*10} FINAL EXPERT SYSTEM REPORT – PATIENT #{patient_number} {'='*10}\n")

        categorized_results = {
            "Genetic": [],
            "Immune": [],
            "Therapy": [],
            "Survival": [],
            "ML": [],
            "Conflict": [],
            "Exceptional": [],
            "Stable": [],
            "Recommendation": [],
        }

        label_map = {
            "Genetic": " Genetic Insights:",
            "Immune": " Immune Profile:",
            "Therapy": " Therapy Evaluation:",
            "Survival": " Survival Insights:",
            "ML": " ML Predictions:",
            "Conflict": " Conflicts:",
            "Exceptional": " Exceptional Cases:",
            "Stable": " Stable Findings:",
            "Recommendation": " Recommendations:",
        }

        def format_certainty(cf):
            percent = round(cf * 100)
            return f"Certainty: {percent}%"

        seen_messages = set()
        results = [f for f in self.facts.values() if isinstance(f, Result)]

        for result in sorted(results, key=lambda x: x.get("certainty_factor", 0), reverse=True):
            r_type = result.get("type")
            msg = result.get("message")
            cf = result.get("certainty_factor", None)

            if r_type in categorized_results and msg not in seen_messages:
                seen_messages.add(msg)
                if cf is not None:
                    categorized_results[r_type].append(f"{msg} [{format_certainty(cf)}]")
                else:
                    categorized_results[r_type].append(msg)

        for category, messages in categorized_results.items():
            if messages:
                print(label_map.get(category, category))
                for msg in messages:
                    print("-", msg)
                print()

        if not any(categorized_results.values()):
            print(" No specific rule matched — recommend full clinical review.\n")

        print(f"{'='*10} End of Report {'='*10}\n")
