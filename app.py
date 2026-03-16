import streamlit as st
import pandas as pd
import joblib


def recommend_treatment(row):
    response = row["Predicted_Response"]
    kras = str(row["KRAS_status"]).lower()
    braf = str(row["BRAF_status"]).lower()
    msi = str(row["MSI_status"]).lower()
    treatment = str(row["treatments.treatment_type"]).lower()

    if "msi" in msi and ("unstable" in msi or "high" in msi):
        return "Immunotherapy (e.g. Pembrolizumab or Nivolumab)"
    if "mut" in braf or "v600e" in braf:
        return "BRAF-targeted combo: Vemurafenib + Cetuximab + Irinotecan"
    if "mut" in kras or "positive" in kras or "yes" in kras:
        return "Avoid EGFR inhibitors; Use FOLFIRI + Bevacizumab"
    if response in [0, 1]:
        return "Escalate treatment: FOLFIRI + Bevacizumab"
    if response in [2, 3]:
        return "Continue FOLFOX or switch to maintenance"
    if response in [4, 5] and "none" not in treatment:
        return "Consider Clinical Trial or Maintenance"
    return "Refer to MDT (Multidisciplinary Team) review"


def generate_ai_explanation(row):
    treatment = row["Recommended_Treatment"]
    response = row["Predicted_Response"]
    notes = []
    if "Immunotherapy" in treatment:
        notes.append("MSI-High status suggests benefit from immunotherapy.")
    elif "BRAF-targeted" in treatment:
        notes.append("BRAF V600E mutation detected; targeted therapy is advised.")
    elif "Avoid EGFR" in treatment:
        notes.append("KRAS mutation found; EGFR inhibitors not effective.")
    elif "Escalate" in treatment:
        notes.append("Low predicted response; escalation recommended.")
    elif "Continue FOLFOX" in treatment:
        notes.append("Good response; continue current protocol.")
    elif "Clinical Trial" in treatment:
        notes.append("Prior treatment + response suggest clinical trial.")
    else:
        notes.append("Complex case. Recommend MDT discussion.")
    return " ".join(notes)





# تحميل النموذج المدرب
model = joblib.load("trained_response_predictor.pkl")

# تصميم الواجهة
st.set_page_config(page_title="Colon Cancer AI Recommender", page_icon="💊", layout="centered")

st.markdown("""
    <style>
    .main {
        background-color: #f0f8ff;
    }
    h1 {
        color: #2a4d69;
        text-align: center;
    }
    </style>
""", unsafe_allow_html=True)

st.title("💊 AI Medical Recommender for Colon Cancer")
st.subheader("🔍 أدخل بيانات المريض:")

# عناصر الإدخال
msi = st.selectbox("MSI status", ["MSI-High", "MSI-Low", "Stable"])
braf = st.selectbox("BRAF status", ["Wild-Type", "V600E", "Mutated"])
kras = st.selectbox("KRAS status", ["Wild-Type", "Mutated", "Positive"])
treatment_type = st.selectbox("Previous Treatment", ["None", "FOLFOX", "FOLFIRI", "Immunotherapy", "Targeted Therapy"])

# زر التنفيذ
if st.button("✨ احسب التوصية"):
    input_data = pd.DataFrame([{
        "MSI_status": msi,
        "BRAF_status": braf,
        "KRAS_status": kras,
        "treatments.treatment_type": treatment_type
    }])

    # التنبؤ
    pred = model.predict(input_data)[0]
    input_data["Predicted_Response"] = pred

    # التوصية والتفسير
    input_data["Recommended_Treatment"] = input_data.apply(recommend_treatment, axis=1)
    input_data["AI_Explanation"] = input_data.apply(generate_ai_explanation, axis=1)

    # العرض
    st.success(f"📈 التنبؤ بالاستجابة: {pred}")
    st.info(f"💉 التوصية: {input_data['Recommended_Treatment'].iloc[0]}")
    st.markdown(f"🧠 تفسير AI: {input_data['AI_Explanation'].iloc[0]}")
