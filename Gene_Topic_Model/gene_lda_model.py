# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 12:57:50 2018

@author: Kari
"""

# Import libraries
import numpy as np
import pandas as pd
import gensim
from gensim import corpora
from gensim.models import CoherenceModel
import matplotlib.pyplot as plt
import pylab as pl

def model_lda(doc_term_lst, doc_names, num_top, num_passes):
    
    # Creating the term dictionary of our courpus, where every unique term is assigned an index. 
    doc_dict = corpora.Dictionary(doc_term_lst)
    
    # Converting list of documents (corpus) into Document Term Matrix using dictionary prepared above.
    doc_term_matrix = [doc_dict.doc2bow(doc) for doc in doc_term_lst] 
    
    # Creating the object for LDA model using gensim library
    Lda = gensim.models.ldamodel.LdaModel
    
    # Running and Trainign LDA model on the document term matrix.
    ldamodel = Lda(doc_term_matrix, num_topics = num_top, id2word = doc_dict, passes = num_passes)
    
    log_perplex = ldamodel.log_perplexity(doc_term_matrix)
#    perplex = 0
    
    tmp_perplex = ldamodel.bound(doc_term_matrix)

    word_perplex = np.exp2(-tmp_perplex / sum(cnt for document in doc_term_matrix for _, cnt in document))

#    print('Coher')
     
#    coher_model = CoherenceModel(model=ldamodel, texts=doc_term_lst, dictionary=doc_dict, coherence='c_v')
#    print('Coher Model Done')
#    coher = coher_model.get_coherence()
#    coher = 0
#    print('Coher Done')

    # Print out term weights per topic
    topic_terms = ldamodel.print_topics(num_topics = num_top, num_words = len(doc_term_lst[0]))
    
    # Create a topic term matrix containing the weights of each term per topic
    terms_lst = list(doc_dict.token2id.keys())
    topic_term_np = np.zeros([len(topic_terms), len(terms_lst)])
    for j in range(len(topic_terms)):
        entry = topic_terms[j]
        temp_lst = entry[1].split('+')
        for i in range(len(temp_lst)):
            temp_str = temp_lst[i]
            ndx = temp_str.find('*')
            if i == (len(temp_lst) - 1):
                temp_term = temp_str[(ndx+2):-1]
            else:
                temp_term = temp_str[(ndx+2):-2]
                
            temp_weight = temp_str[0:ndx]
            
            term_ndx = terms_lst.index(temp_term)
            topic_term_np[j,term_ndx] = temp_weight
                
    topic_names = ['topic_'+str(x) for x in range(len(topic_terms))]
    topic_term_pd = pd.DataFrame(topic_term_np, columns = terms_lst, index = topic_names)
    
    # Get doc/topic matrix
    corpus_lda_dense  = gensim.matutils.corpus2dense(ldamodel[doc_term_matrix], num_terms = ldamodel.num_topics, 
                                                     num_docs=doc_dict.num_docs)
    
    # Convert the doc/topic matrix to a panda datagframe
    doc_topic_np = np.array(corpus_lda_dense).transpose()
#    doc_names = ['doc_'+str(x) for x in range(doc_topic_np.shape[0])]
    doc_topic_pd = pd.DataFrame(doc_topic_np, columns = topic_names, index = doc_names)
     
    return doc_topic_pd, topic_term_pd, log_perplex, word_perplex, ldamodel, doc_term_lst, doc_dict

# Matrix passed in is assumed to be rows = genes, columns = cell identifiers or UMIs 
# (UMIs are docs, genes are terms)
def convert_matrix(matrix, gene_names, cell_names):
    
    matrix_np = np.array(matrix)
    matrix_np = matrix_np.transpose()
    
    doc_l2d_st = []
    for i in range(len(cell_names)):
        tmp_ndx = matrix_np[i,:] != 0
        
        tmp_lst = []
        for j in range(len(tmp_ndx)):
            if tmp_ndx[j]:
                tmp_lst.append(gene_names[j])
                
        doc_l2d_st.append(tmp_lst)
        
    return doc_l2d_st   
    
def perform_lda(matrix, gene_names, cell_names, topic_vals):
    
    doc_term_lst = convert_matrix(matrix, gene_names, cell_names)
    
    log_perplex_lst = []
    word_perplex_lst = []
    for n in topic_vals:
        
        doc_topic_pd, topic_term_pd, log_perplex, word_perplex, ldamodel, doc_term_lst, doc_dict =  model_lda(doc_term_lst, cell_names, n, 1)
        log_perplex_lst.append(log_perplex)
        word_perplex_lst.append(word_perplex)
    
    return doc_topic_pd, topic_term_pd, log_perplex_lst, word_perplex_lst, ldamodel, doc_term_lst, doc_dict
    
    
file_name = 'C:\\DePaulCoursework\\Research\\Data\\Sample1_Infected_5000UMI_Example.csv'
gene_pd = pd.read_csv(file_name, delimiter = ',', index_col = 0, header = 0)
print(gene_pd.head())

cell_names = list(gene_pd.columns.values)
gene_names = list(gene_pd.index)

gene_np = np.array(gene_pd)

topic_vals = list(np.arange(10,21))
#topic_vals = [1]
doc_topic_pd, topic_term_pd, log_perplex_lst, word_perplex_lst, ldamodel, doc_term_lst, doc_dict = perform_lda(gene_np, gene_names, cell_names, topic_vals)

print('Doc - Topic Proportion DataFrame:')
print(doc_topic_pd.head())
print('\nTopic - Term Weight DataFrame:')
print(topic_term_pd.head())
 
pl.figure()
pl.xlabel("Number of Topics")
pl.ylabel("Log Perplexity")
pl.title("Log Perplexity Vs Number of Topics")
pl.plot(topic_vals, log_perplex_lst)
plt.show()
    
pl.figure()
pl.xlabel("Number of Topics")
pl.ylabel("Per Word Perplexity")
pl.title("Per Word Perplexity Vs Number of Topics")
pl.plot(topic_vals, word_perplex_lst)
plt.show()

    