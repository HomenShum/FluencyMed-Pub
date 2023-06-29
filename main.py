import openai
import re
import tenacity
import time
import os
import pinecone
from langchain.embeddings.openai import OpenAIEmbeddings
from langchain.vectorstores import Pinecone
import streamlit as st
from streamlit_chat import message
from Bio import Entrez

##### Settings ############################################################################################################
### Use streamlit to ask user input for openai.api_key
openai.api_key = st.text_input("Enter your OpenAI API Key", value="", type="password")

if openai.api_key:
    # Pinecone Settings
    os.environ['PINECONE_API_KEY'] = 'd0dc3dbd-8599-4d68-9816-6a0d74a290c5'
    os.environ['PINECONE_ENVIRONMENT'] = 'us-west4-gcp'
    pinecone.init(
        api_key=os.getenv('PINECONE_API_KEY'), 
        environment=os.getenv('PINECONE_ENVIRONMENT')
    )
    index_name = 'langchainpdfchat'
    embeddings = OpenAIEmbeddings()
    docsearch = Pinecone.from_existing_index(index_name, embeddings)



    ##### Functions ############################################################################################################
    @tenacity.retry(
        stop=tenacity.stop_after_delay(30),
        wait=tenacity.wait_exponential(multiplier=1, min=1, max=30),
        retry=tenacity.retry_if_exception_type(openai.error.APIError),
        reraise=True,
    )
    def gpt_completion(
        prompt,
        engine="gpt-3.5-turbo-0613",
        temp=0,  # set at 0 to ensure consistent completion, increase accuracy along with UUID
        top_p=1.0,
        tokens=1000,
        freq_pen=0.25,
        pres_pen=0.0,
        stop=["<<END>>"],
    ):
        prompt = prompt.encode(encoding="ASCII", errors="ignore").decode()
        response = openai.ChatCompletion.create(
            model=engine,
            messages=[
                {"role": "system", "content":   
                        """You are an AI assistant specialized in biomedical topics. You are provided with a text description from a patient's screening notes. Analyze the patient's notes and ask follow up question. Here are your instructions:

                        - Highlight medical advice or diagnostic information from sample dialogues dataset. 

                        - Ensure the output is in markdown bullet point format for clarity.

                        - Encourage the user to consult a healthcare professional for advice."""},
                {"role": "user", "content": prompt},
            ],
            max_tokens=tokens,
            temperature=temp,
            top_p=top_p,
            frequency_penalty=freq_pen,
            presence_penalty=pres_pen,
            stop=stop,
        )
        text = response["choices"][0]["message"]["content"].strip()
        text = re.sub("\s+", " ", text)
        return text

    def analyze_patient_note(patient_note):
        start_time = time.time()

        prompt = """- Ask the patient follow up diagnosis question using the patient notes and Clinical Classifications Software Refined (CCSR) categories. 
                    
                - Highlight medical advice or diagnostic information from ~370,000 sample patient physician dialogues dataset. 

                - The output should be detailed and gives physicians a top down view of the medical situation quickly. 
                
                - Encourage the user to consult a healthcare professional for advice.

                \n
                """

        input_var_1 = "Patient Notes: \n"

        input_patient_note_analysis = patient_note_analysis(patient_note)

        input_var_2 = "Clinical Classifications Software Refined (CCSR) categories listed below: \n"

        prompt = prompt + input_var_1 + input_patient_note_analysis + input_var_2

        ccsr_categories_list = search_index(input_patient_note_analysis, ccsr_df_feather['Embeddings'].tolist())

        ccsr_categories_list_list = []

        for i, category in enumerate(ccsr_categories_list, start=1):
            content = category['content']
            prompt += f"{i}. {content}\n"
            ccsr_categories_list_list.append(content)

        search = docsearch.similarity_search(prompt)
        input_dialogues_data = f"sample physician diagnosis dataset: \n{search[0].page_content}\n"
        prompt = prompt + input_dialogues_data +    """            
                                                    Output in markdown or list format:
                                                    **Medical advice or diagnostic information based on similar real diagnostics:**
                                                    <<Output Here>>
                                                    """

        result = gpt_completion(prompt)

        end = time.time()
        print(f"Runtime of the program is {end - start_time}")

        return result, input_patient_note_analysis, ccsr_categories_list_list

    def get_patient_notes():
        input_text = st.text_area("Patient Notes",
                                "NEED HELP. PROLONGED FEVER WITH NO OTHER SYMPTOMS. For history : I have a 1 year old son, 10kg , who started day care a few months ago. We've been consistently sick with a bunch of stuff since. Gastro, flu colds etc. Last weekend he started vomiting and having a fever. We took him to the ER because he started shaking pretty badly in his sleep, almost convulsing. We couldnt tell what. He was out of it and not really responsive when we picked him up and was just looking off into the distance. Anyway. They recorded a fever of 40 C/ 104 F. Took us in and gave nurofen, and did a piss test and found nothing. Discharged us when the fever went down. He stopped vomiting the next day. A few days go past and he's better but still has fevers randomly. One night I notice that his left foot is sore. And he doesn't want to use it. I also noticed he hasn't really been walking much. Another day or 2 goes by and his foot is better. But now his left shoulder appears to be injured. He won't lift his arm above shoulder height and again doesn't want to walk because he can't use his left arm to balance. A doctor examined him on Wednesday. Says his arm isn't broken or dislocated. It's just a random injury and doesn't comment on the foot. The doctor says teething can cause low grade fevers and sends us home. Arm gets better the next day and is fine by Friday. Today, 8 days after the initial ER visit. He's having 100F or 38 C fevers randomly throughout the day. He has zero other symptoms. He's walking again. His arm is fine again. His foot is normal. It's been probably 6 or 7 days since he's had any sick symptoms like vomiting or runny nose. I'm going to the doctor tomorrow with him. But doctors are generally useless in my town. Need help or suggestions. Is this normal?" ,
                                    key = 'patient_notes',
                                    height = 300)
        return input_text

    def search(query):
        Entrez.email = 'hshum2018@gmail.com'
        handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax='1',
                                retmode='xml',
                                term=query)
        results = Entrez.read(handle)
        return results

    def fetch_details(id_list):
        ids = ','.join(id_list)
        Entrez.email = 'hshum2018@gmail.com'
        handle = Entrez.efetch(db='pubmed',
                            retmode='xml',
                            id=ids)
        results = Entrez.read(handle)
        return results

    def get_abstract(paper):
        abstract = ''
        if 'Abstract' in paper['MedlineCitation']['Article']:
            abstract = paper['MedlineCitation']['Article']['Abstract']['AbstractText']
            if isinstance(abstract, list):
                abstract = ' '.join(abstract)
        return abstract

    def format_text(text):
        formatted_text = "\n".join(f"- {item.strip()}" for item in text.split("-")[1:])
        if "diagnostic information" in text.lower():
            return "Medical advice or diagnostic information based on similar real diagnostics:\n" + formatted_text
        else:
            return formatted_text


    ##### Main #####################################
    st.sidebar.header("Recommended PubMed Articles")

    st.title("Fluency Med: Analyze Patient Notes and Ask Follow Up Diagnosis Question")
    st.write("This is a chatbot that can answer questions about medical advice or diagnostic information from sample dialogues dataset. It is trained on the GPT-3.5-Turbo engine and is specialized in biomedical topics. It is provided with a text description from a patient's screening notes. Analyze the patient's notes and ask follow up question. Here are your instructions:")
    st.write("- Highlight medical advice or diagnostic information from sample dialogues dataset.")
    st.write("- Ensure the output is in markdown bullet point format for clarity.")
    st.write("- Encourage the user to consult a healthcare professional for advice.")

    st.write("Please enter your patient's screening notes below:")

    ### storing the chat
    if 'generated' not in st.session_state:
        st.session_state['generated'] = []

    if 'past' not in st.session_state:
        st.session_state['past'] = []

    patient_note = get_patient_notes()

    # Add reset chat button
    if st.button('Reset Chat'):
        st.session_state['generated'] = []
        st.session_state['past'] = []

    # Analyze patient notes button
    if st.button('Analyze Patient Notes'):  

        if openai.api_key and patient_note:  # Only run if both API key and patient notes are provided

            from feature2_ccsr_categorization import search_index, ccsr_df_feather
            from feature1_clinical_note_summarization import patient_note_analysis
            
            # Run GPT model to analyze patient's note
            result, input_patient_note_analysis, ccsr_categories_list_list = analyze_patient_note(patient_note)

            # Format the result and CCSR categories as a bulleted list
            formatted_result = format_text(result)
            formatted_ccsr_categories = "CCSR Categories Identified:\n" + "\n".join(f"- {item}" for item in ccsr_categories_list_list)

            # Store the output
            st.session_state.past.append(patient_note) 
            st.session_state.generated.append(formatted_result)  # Store the formatted result
            st.session_state.generated.append(input_patient_note_analysis)
            st.session_state.generated.append(formatted_ccsr_categories)  # Store the formatted categories

            if st.session_state['generated']:
                for i in range(len(st.session_state['past'])-1, -1, -1):
                    message(st.session_state['past'][i], is_user=True, key=str(i) + '_user')
                    message(st.session_state['generated'][3*i], key=str(3*i))  # Display the formatted result
                    message(st.session_state['generated'][3*i+1], key=str(3*i+1))
                    message(st.session_state['generated'][3*i+2], key=str(3*i+2))  # Display the formatted categories

            # Perform PubMed search for each category
            for category in ccsr_categories_list_list:
                results = search(category)
                id_list = results['IdList']
                papers = fetch_details(id_list)
                for i, paper in enumerate(papers['PubmedArticle']):
                    title = paper['MedlineCitation']['Article']['ArticleTitle']
                    author_list = paper['MedlineCitation']['Article']['AuthorList']
                    authors = ', '.join([author.get('LastName', '') for author in author_list])
                    abstract = get_abstract(paper)
                    
                    # Show article in sidebar
                    st.sidebar.subheader(f'Title: {title}')
                    st.sidebar.write(f'Authors: {authors}')
                    st.sidebar.write(f'Abstract: {abstract}')
