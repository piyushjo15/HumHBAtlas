## this script is for Word to vec comparing the graf network obtained from
## normal glial cells to glial cels in PA tumors as an example
## similar script can be run for DMG samples
import os
import sys
import networkx as nx
from gensim.models import Word2Vec
from scipy.spatial.distance import cosine
from scipy.linalg import orthogonal_procrustes
import pandas as pd
import numpy as np
import pickle
from scipy.sparse import load_npz
from scipy.stats import pearsonr
from sklearn.metrics.pairwise import cosine_similarity
import pecanpy as pc
## set seed for reproducibility
os.environ['PYTHONHASHSEED'] = '42'

## ALL samples can be comapred together

proj_DIR="/home//RECORDR/"
os.chdir(proj_DIR)
print("Performing Word2Vec analysis for PA samples ")

## 1.1 load the graphs and generate node2vec using pecanpy
## pecanpy function
def run_node2vec(edg_path, p=1, q=1, num_walks=200, walk_length=80, workers=10, seed=42):
    n2v = pc.pecanpy.DenseOTF(p=p, q=q, workers=workers, random_state=seed)
    n2v.read_edg(edg_path, weighted=True, directed=False, delimiter='\t')
    walks = n2v.simulate_walks(num_walks=num_walks, walk_length=walk_length)
    return walks

##1.2 Train Word2Vec model

def train_w2v(walks, vector_size=64, window=10, min_count=5, sg=1, negative=5, epochs=10, workers=10, seed=42):
    return Word2Vec(
        sentences=walks, vector_size=vector_size, window=window,
        min_count=min_count, sg=sg, negative=negative, alpha=0.025,
        epochs=epochs, workers=workers, seed=seed
    )

# Normal
## there reference are created for normal
partx='a' # 'b' or 'c'

walks_nor= run_node2vec(f"Normalv2/gene_network_nor_{partx}.edg")
w2v_nor= train_w2v(walks_nor)

# Tumors
walks_tum_PA01= run_node2vec("PA01/gene_network_PA01.edg")
w2v_tum_PA01= train_w2v(walks_tum_PA01)

walks_tum_PA02= run_node2vec("PA02/gene_network_PA02.edg")
w2v_tum_PA02= train_w2v(walks_tum_PA02)

walks_tum_PA03= run_node2vec("PA03/gene_network_PA03.edg")
w2v_tum_PA03= train_w2v(walks_tum_PA03)



##align W2V using Amir's script

#smart procrustes alignment functions - taken from - https://gist.github.com/zhicongchen/9e23d5c3f1e5b1293b16133485cd17d8 - this is freely available - for gensim model alignment (gensim version >= 4) \
#based on https://arxiv.org/pdf/1605.09096.pdf
def smart_procrustes_align_gensim(base_embed, other_embed, words=None):
    """
    Original script: https://gist.github.com/quadrismegistus/09a93e219a6ffc4f216fb85235535faf
    Procrustes align two gensim word2vec models (to allow for comparison between same word across models).
    Code ported from HistWords <https://github.com/williamleif/histwords> by William Hamilton <wleif@stanford.edu>.
        
    First, intersect the vocabularies (see `intersection_align_gensim` documentation).
    Then do the alignment on the other_embed model.
    Replace the other_embed model's syn0 and syn0norm numpy matrices with the aligned version.
    Return other_embed.

    If `words` is set, intersect the two models' vocabulary with the vocabulary in words (see `intersection_align_gensim` documentation).
    """

    # patch by Richard So [https://twitter.com/richardjeanso) (thanks!) to update this code for new version of gensim
    # base_embed.init_sims(replace=True)
    # other_embed.init_sims(replace=True)

    # make sure vocabulary and indices are aligned
    in_base_embed, in_other_embed = intersection_align_gensim(base_embed, other_embed, words=words)

    # get the (normalized) embedding matrices
    base_vecs = in_base_embed.wv.get_normed_vectors()
    other_vecs = in_other_embed.wv.get_normed_vectors()

    # just a matrix dot product with numpy
    m = other_vecs.T.dot(base_vecs) 
    # SVD method from numpy
    u, _, v = np.linalg.svd(m)
    # another matrix operation
    ortho = u.dot(v) 
    # Replace original array with modified one, i.e. multiplying the embedding matrix by "ortho"
    other_embed.wv.vectors = (other_embed.wv.vectors).dot(ortho)    
    
    return other_embed

def intersection_align_gensim(m1, m2, words=None):
    """
    Intersect two gensim word2vec models, m1 and m2.
    Only the shared vocabulary between them is kept.
    If 'words' is set (as list or set), then the vocabulary is intersected with this list as well.
    Indices are re-organized from 0..N in order of descending frequency (=sum of counts from both m1 and m2).
    These indices correspond to the new syn0 and syn0norm objects in both gensim models:
        -- so that Row 0 of m1.syn0 will be for the same word as Row 0 of m2.syn0
        -- you can find the index of any word on the .index2word list: model.index2word.index(word) => 2
    The .vocab dictionary is also updated for each model, preserving the count but updating the index.
    """

    # Get the vocab for each model
    vocab_m1 = set(m1.wv.index_to_key)
    vocab_m2 = set(m2.wv.index_to_key)

    # Find the common vocabulary
    common_vocab = vocab_m1 & vocab_m2
    if words: common_vocab &= set(words)

    # If no alignment necessary because vocab is identical...
    if not vocab_m1 - common_vocab and not vocab_m2 - common_vocab:
        return (m1,m2)

    # Otherwise sort by frequency (summed for both)
    common_vocab = list(common_vocab)
    common_vocab.sort(key=lambda w: m1.wv.get_vecattr(w, "count") + m2.wv.get_vecattr(w, "count"), reverse=True)
    # print(len(common_vocab))

    # Then for each model...
    for m in [m1, m2]:
        # Replace old syn0norm array with new one (with common vocab)
        indices = [m.wv.key_to_index[w] for w in common_vocab]
        old_arr = m.wv.vectors
        new_arr = np.array([old_arr[index] for index in indices])
        m.wv.vectors = new_arr

        # Replace old vocab dictionary with new one (with common vocab)
        # and old index2word with new one
        new_key_to_index = {}
        new_index_to_key = []
        for new_index, key in enumerate(common_vocab):
            new_key_to_index[key] = new_index
            new_index_to_key.append(key)
        m.wv.key_to_index = new_key_to_index
        m.wv.index_to_key = new_index_to_key
        
        print(len(m.wv.key_to_index), len(m.wv.vectors))
        
    return (m1,m2)


#cosines for timecourse - orthogonal procrustes
def get_cosines(reference_model, 
                reference_vocab: list, 
                name_list: list, 
                reference_name: str, 
                write_path_ref_model: str,
                *comparison_models):
    '''
    Function to get cosine similarity between different gensim models. takes a reference model and aligns
    them. the model being compared is aligned to the reference model. 
    
    Arguments:
    --------
    reference_model
        gensim model -  model that all other models will be aligned to (in a pairwise manner)
    reference_vocab
        List -  reference vocabulary, should be vocabulary from reference model.
    comparison_models
        gensim models - any number of models can be passed. takes in gensim models and these
        are the models that will be aligned to reference model so that vectors can be compared

    Returns:
    -------
    returns a dataframe of all pairwise comparisons with the cosine similarity of gene 'x' in reference model
    relative to that same gene in comparison model. Only common vocabulary is used when aligning

    Alignment happens using orthogonal procrustes. Reference implementation(freely available) - https://gist.github.com/zhicongchen/9e23d5c3f1e5b1293b16133485cd17d8
    based on paper -  https://arxiv.org/pdf/1605.09096.pdf
    '''
    df = pd.DataFrame(index=[reference_name])

    if not os.path.exists(write_path_ref_model):
        os.makedirs(write_path_ref_model)
    reference_model.save(os.path.join(write_path_ref_model, f'OligoAstroMicroglia_{partx}.model'))

    for name_, w2vmodel in zip(name_list, comparison_models):
        print(f'aligning to {name_}')
        reference = Word2Vec.load(os.path.join(write_path_ref_model, f'OligoAstroMicroglia_{partx}.model'))
        print(f'length of reference vocab: {len(reference.wv.index_to_key)}')

        smart_procrustes_align_gensim(reference, w2vmodel)
        
        ##save the reference and rotated model
        common_genes = [g for g in reference.wv.index_to_key if g in w2vmodel.wv.key_to_index]
        ##obtain embedding matrix
        E_ref = np.vstack([reference.wv[g] for g in common_genes]).astype(np.float32)
        E_tum = np.vstack([w2vmodel.wv[g]  for g in common_genes]).astype(np.float32)

        np.savez(
          os.path.join(os.fspath(name_), f"embeddings_aligned_{partx}.npz"),
          genes=np.array(common_genes),
          ref=E_ref,
          tum=E_tum
          )
        
        cos_dict = {}
        for j in reference_vocab:
            try:
                cos_dict[j] = 1 - cosine(reference.wv[j], w2vmodel.wv[j])# 1 - cosine for cosine similarity: scipy does cosine distance
            except KeyError:
                cos_dict[j] = -1.1 ##gene absent in tumor

        tempdf = pd.DataFrame([cos_dict])
        tempdf.index = [name_]

        df = pd.concat([df, tempdf])
        df.iloc[:1] = 1.0
    try:
        os.remove(os.path.join(write_path_ref_model, f'OligoAstroMicroglia_{partx}.model'))
        shutil.rmtree(write_path_ref_model)
    except Exception as e:
        print(f'Error during cleanup: {e}')
        print('Reference model and folder not deleted, please delete manually')

    return df



out_path=proj_DIR+"Normalv2/"
genes_nor=set(w2v_nor.wv.key_to_index)
df_similarity=get_cosines(
    w2v_nor,
    genes_nor,
    ["PA01","PA02", "PA03"],
    "OligoAstroMicroglia",
    out_path,
    w2v_tum_PA01,
    w2v_tum_PA02,
    w2v_tum_PA03
)


##this output is a data.frame with cosine values for all the genes present in 
## reference. if the gene is absent in tumor it get -1.1 score in tumor
## genes present in tumor are not counted, I should add this letter for each tumor
with open(f"PA_cosines_{partx}.pkl", "wb") as f:
 pickle.dump(df_similarity, f)

exit
