# FluencyMed-Pub
Increasing doctor capacity and patient care quality. We enrich physician practices and enhance patient telehealth experience. Fluency Med is creating a healthier, more efficient future in healthcare.

The code is private for development purposes, for collaboration, please reach out to me at hshum2018@gmail.com

Inspiration
When my family needs to find a doctor for a random symptom that is in dire need to get treatment, I frequently found ourselves turmoiling over the time wasted on booking appointments between each and every part of the medical department. Sometimes, the symptom magically disappears, and we no longer need to attend the appointment. Then sometimes the symptoms would come back, and my family would need immediate medical attention. However, due to our reluctance to face the friction created by the long dragging appointment system for hospitals, we would have to start from ground zero in order to get a primary physician, then get a referral after a series of test and questionnaires, and then not get the right doctor, and then ultimately loop around for at least a few months until we finally find the right department, hear about the right terminologies used to describe our symptoms and issues, and get the medications needed for treatment.

This knowledge gap between patient and physician is too wide and the appointment process to get to know what the heck is wrong is us takes too long. If the symptom already subsided/went to hibernation during the diagnosis process, then it would be left neglected due to the need to take care of other priorities in this increasingly complex and busy society that we are living in.

Therefore, the call to action led me to create this experimental MVP to help a potential patient in dire need to figure out what is going on with them and create a summary that is concise and easy to communicate to the physician using 370K patient-physician dialogues and the 35 million articles publicized on PubMed website.

What it does
Takes in patient's description of everything that is wrong
Analysis and output summary that is easy to read in order to find the most important information: such as the patient information, symptoms and conditions, concerns, medical history and investigations, and recommended medical department to contact.
Use the output to generate Clinical Classifications Software Refined (CCSR) categories, which is then used to search for relevant articles from PubMed that can supplement patient and physicians on knowledge needed to understand the problem that they encounter.
The output is also used to perform a similarity search on a pinecone database that has 370K total patient and physician diagnosis dialogues. The most similar dialogue is then used to help understand the situation and quote out any medical diagnosis that would accelerate the identification of the patient's symptoms.
How we built it
GPT 3.5 Turbo is used to generate the analysis summary for all aspect of the application
Used llama_index to upsert the embeddings of the 370K dialogues onto Pinecone Vector Database
Used langchain to perform the similarity search using the patient note analysis on the Pinecone database
Used Streamlit to host the user interface.
Challenges we ran into
The existing database would need further investment and development in order to scale larger and be updated with the more latest information available for medical diagnosis support. Existing solution to pull full text articles from PubMed is poor, and so I was not able to extract key information from the recommended full text articles to deduce the research time needed to get accessible medical domain knowledge tailored to patient need, yet.

Accomplishments that we're proud of
Proud of the application's capability to understand and retrieve the most relevant information using real world case scenario to help patients with their urgent needs and assist in communication.

What we learned
It would not be feasible to use GPT to perform diagnosis directly as it can push away the physicians, especially independent practices, since it would take over their jobs. Instead, a supportive function that can enhance patient-physician communication and knowledge is more advanced in the need of the current society.

What's next for Fluency Med x Semantic Scholar
To further develop the generative AI use cases in the medical field and better assist in the physician side of things to save their time and allow them ability to provide better cares to the patients!

## Future Roadmap:

Business Model:

The business model for the Fluency Med services is a dual strategy involving both a paid subscription model for patients and a business-to-business (B2B) model for healthcare providers. Services include billing, documentation, administrative work, appointment processes, and diagnostic assistance via semantic search on knowledge database. We focus on selling the service to independent providers, associations, and healthcare software providers, and we help affluent individuals with their specific needs on getting to know their medical needs better while making their process more efficient.

Proof of Market Fit:

Fluency Med services addresses several critical challenges in the healthcare industry, demonstrating a strong market fit. The service helps to reduce physician burnout by streamlining administrative tasks, a common issue in healthcare. It improves patient outcomes by offering diagnostic support and enhancing the online care user experience. By improving physician satisfaction and attracting more medical students to primary care, the service demonstrates that it meets both industry and user needs. Further evidence of market fit comes from its affordability, which appeals to cost-conscious consumers, and the comprehensive nature of its solutions, which fulfill the desire for a one-stop healthcare service.

Traction:

Traction for the Fluency Med services is demonstrated by the numerous activities and initiatives already in progress. Interviews between patients and providers, discussions with pilot customers, and online marketing efforts have generated interest and potential customers. The pilot group, the Christian healthcare group, consisting of 180 physicians serving an estimated patient base of 6,000, showcases the service's operational feasibility. Furthermore, collaborations with healthcare and incubator organizations to help in securing alpha customers, letters of intent, and fund-raising, showing that the service is gaining momentum in the market.

Starting with the patient interface, we can go through the key features.

Patient

GPT Summary and analysis
Create CCSR categories
Leading Questions supported by 382K patient-doctor conversations dataset that are real plus generated using OpenAI’s GPT LLM.
The datasets are as follows: en_medical_dialog (200K), GenMedGPT-5k, HealthCareMagic-100k, iCliniq
LLAMA-INDEX and PINECONE vector database for index search
End of patient screening session notes full-analysis using all dialogues and patient information, treatment recommendation, medication recommendation, and preliminary diagnosis is supported by Clinical Classifications Software Refined (CCSR) categorization as well as real patient & physician diagnostic notes.
Patient Physician

OpenAI transcription and secured live call feature
GPT summary and analysis for live transcription
Physician

Supplemental Medical Research from PubMed 35million biomedical research articles.
Final diagnosis note signed off by physician is used to find top billing code suggestions with evidence provided to enhance administrative process flow
Appointment process management for physician to help bump off none showing patients quickly and append patients that want faster care.
CosmaNeura Pitch Deck 16 June 23.pptx.pdf

Demo Physician Interface:

PhysicianSideCalHacksCosmas

Demo Patient Interface:

CalHackCosmas

GPT Prompts: messages=[ { "role": "system", "content": """You are an AI assistant specialized in biomedical topics. You are provided with a text description from a patient's screening notes. Analyze the patient's notes and provide information useful for Clinical Classifications Software Refined (CCSR) categorization. Here are your instructions:

            - Highlight conditions, symptoms, medical history, and any other information that can be mapped to specific CCSR categories.

            - Keep in mind that the CCSR is used for grouping a large number of diseases into manageable categories for statistical analysis and reporting. 

            - Ensure the conversation includes information that can guide the mapping to CCSR categories, such as the type of disease, cause, location in the body, and patient's age and sex. 

            - Avoid providing medical advice or diagnostic information. 

            - Ensure the output is in markdown bullet point format for clarity.

            - Encourage the user to consult a healthcare professional for advice."""
        },
        {"role": "user", "content": "Can you analyze patient 123's notes and provide useful information for CCSR categorization? Output in markdown bullet points."}
    ],

    messages=[
        {"role": "system", "content":   """You are an AI assistant specialized in biomedical topics. You are provided with a text description from a patient's screening notes. Analyze the patient's notes and ask follow up question. Here are your instructions:

                - Highlight conditions, symptoms, medical history, and any other information that can be mapped to specific CCSR categories.

                - Keep in mind that the CCSR is used for grouping a large number of diseases into manageable categories for statistical analysis and reporting. 

                - Ensure the conversation includes information that can guide the mapping to CCSR categories, such as the type of disease, cause, location in the body, and patient's age and sex. 

                - Highlight medical advice or diagnostic information quoted and summarized from the given information. 

                - Ensure the output is in markdown bullet point format for clarity.

                - Encourage the user to consult a healthcare professional for advice."""},

        {"role": "user", "content": prompt},
    ],
prompt = """Provide an end session screening summary and provide the most recommended medical department for the patient to contact, provide the patient's preliminary medical diagnosis, emphasize on "why" for the diagnosis.

        - Highlight or quote any medication and treatment recommendation based on the patient note and sample dialogues dataset. 
        
        - The output should be appropriate for the patient's age. 
        
        - Ensure the output is in markdown bullet point format for clarity.

        - Encourage the user to consult a healthcare professional for advice.
        \n
        """
input_var_1 = "Patient Notes:"

input_end_session_note = """

                    "input": "Hello doctor,I had mumps five months ago and after that, I started to have an infection in my left testes. It was swollen and now it has shrunk to almost half the size of the other one. As I am sexually active, I feel a pain in each of the vas deferens after sex. If I do not have sex for days, they become sensitive. I was treated with Ceftum 500 mg, the first time I had an infection. Now my question is, is there any chance that the infection is still in my body? And, do I need to get examined for it? For the time being, please suggest some precautionary antibiotics for my relief.",

                    \n
                    """
input_var_2 = "Clinical Classifications Software Refined (CCSR) categories listed below: \n"

prompt = prompt + input_var_1 + input_end_session_note + input_var_2

for i, category in enumerate(ccsr_categories_list, start=1): content = category['content'] prompt += f"{i}. {content}\n"

input_dialogues_data = "sample dialogues dataset: \n"

input_var_3 = res

prompt = prompt + input_dialogues_data + input_var_3

##Sample Real P&P Dialogue## """ "patient": "Hello doctor,My son aged 13 years has got his left eye -11 and right with -8. I am now very worried what good precautions should I take with him. Can you suggest me what should I do?", "physician": "Hello, Welcome to Chat Doctor forum. I have gone through what you have mentioned and understand your anxiety about your son\s eyesight. - 11D and - 8D are high myopic powers and the chances are that the power may progress as your son is still young. My advice to you is: 1. Avoid reading or holding objects very near to the eye. This will cause myopia to progress faster and will worsen the condition. 2. Every six monthly visits or follow up visits as per doctors advise. Since the myopia is of a high degree, there will be excessive stretching and degeneration of the structures of the eyeball. This can cause problems such as retinal detachment. Regular follow-ups help us in picking up changes faster and in preventing such complications as retinal detachment. 3. Ask your son to regularly self-check his vision. This is done by closing each eye one at a time and observing how a particular object or letter looks. In case he is able to pick up any change or ChatDoctor. 4. Urgently visit your doctor in case of a sudden ", """
