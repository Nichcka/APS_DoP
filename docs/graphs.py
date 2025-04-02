import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from datetime import datetime

def load_data(input_file):
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Файл {input_file} не найден!")
    df = pd.read_csv(input_file, sep='\t')
    data_type = 'protein' if 'Per_Sim' in df.columns else 'dna'
    print(f"Тип данных: {data_type.upper()}")
    return df, data_type

def create_matrices(df, data_type):
    all_ids = sorted(set(df['Seq_1']).union(set(df['Seq_2'])))
    if data_type == 'dna':
        identity_matrix = pd.DataFrame(100, index=all_ids, columns=all_ids)
        for _, row in df.iterrows():
            id1, id2 = row['Seq_1'], row['Seq_2']
            identity_matrix.loc[id1, id2] = identity_matrix.loc[id2, id1] = row['Per_Id']
        return identity_matrix, None
    else:
        identity_matrix = pd.DataFrame(100, index=all_ids, columns=all_ids)
        similarity_matrix = pd.DataFrame(100, index=all_ids, columns=all_ids)
        for _, row in df.iterrows():
            id1, id2 = row['Seq_1'], row['Seq_2']
            identity_matrix.loc[id1, id2] = row['Per_Id']
            similarity_matrix.loc[id1, id2] = row['Per_Sim']
        return identity_matrix, similarity_matrix

def plot_heatmaps(id_matrix, sim_matrix, data_type):
    os.makedirs("results", exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # График идентичности
    matrix_size = id_matrix.shape[0]
    figsize = (max(10, matrix_size * 0.8), max(8, matrix_size * 0.6))  # Adjust factors
    plt.figure(figsize=figsize)
    sns.heatmap(id_matrix, annot=True, fmt=".1f", cmap="YlOrRd", vmin=0, vmax=100)
    plt.title("Pairwise Identity (%)")
    identity_file = f"results/identity_{timestamp}.png"
    plt.savefig(identity_file, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Сохранен: {identity_file}")
    
    # График сходства (только для белков)
    if data_type == 'protein':
        matrix_size = id_matrix.shape[0]
        figsize = (max(10, matrix_size * 0.8), max(8, matrix_size * 0.6))  # Adjust factors
        plt.figure(figsize=figsize)
        sns.heatmap(sim_matrix, annot=True, fmt=".1f", cmap="YlGnBu", vmin=0, vmax=100)
        plt.title("Pairwise Similarity (%)")
        similarity_file = f"results/similarity_{timestamp}.png"
        plt.savefig(similarity_file, bbox_inches='tight', dpi=300)
        plt.close()
        print(f"Сохранен: {similarity_file}")

