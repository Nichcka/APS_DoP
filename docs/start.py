import argparse

# Крч сюда вставляем основной:
# def process_alignment_DNA(input_file, output_file):
# def process_alignment_protein(input_file, output_file, aminoacid):

def file_exists(filepath):
    try:
        with open(filepath, 'r'):
            return True
    except FileNotFoundError:
        return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-m', '--mode')
    parser.add_argument('-o', '--output')
    parser.add_argument('-amk', '--aminoacid')
    parser.add_argument('--dpi', type=int, default=300, help='Качество изображений (по умолчанию: 300)')
    args = parser.parse_args()

    if args.mode == 'DNA':
        process_alignment_DNA(args.input, args.output)
    elif args.mode == 'protein':
        process_alignment_protein(args.input, args.output, args.aminoacid)
    
    df, data_type = load_data(args.-i)
    id_matrix, sim_matrix = create_matrices(df, data_type)
    plot_heatmaps(id_matrix, sim_matrix, data_type)
