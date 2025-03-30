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
    args = parser.parse_args()

    if args.mode == 'DNA':
        process_alignment_DNA(args.input, args.output)
    elif args.mode == 'protein':
        process_alignment_protein(args.input, args.output, args.aminoacid)
