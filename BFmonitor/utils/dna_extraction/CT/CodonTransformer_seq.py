#modules for CodonTransformer
import torch
import argparse
from transformers import AutoTokenizer, BigBirdForMaskedLM
from CodonTransformer.CodonPrediction import predict_dna_sequence


def CodonTransformer_process(protein, organism):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    #Argument parsing
    parser=argparse.ArgumentParser()
    parser.add_argument('--protein', '-p' , required=True, help="Sequence of the protein whose want to be reverse translate")
    parser.add_argument('--organism', '-o', default='Escherichia coli general', help="Organism in which the protein is going to be expressed")

    args=parser.parse_args()


    # Load model and tokenizer
    tokenizer = AutoTokenizer.from_pretrained("adibvafa/CodonTransformer")
    model = BigBirdForMaskedLM.from_pretrained("adibvafa/CodonTransformer").to(device)


    # Set your input data
    protein = args.protein
    organism = args.organism


    # Predict with CodonTransformer
    output = predict_dna_sequence(
        protein=protein,
        organism=organism,
        device=device,
        tokenizer=tokenizer,
        model=model,
        attention_type="original_full",
        deterministic=True
    )

    dna_seq=output.predicted_dna

    return dna_seq

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--protein', '-p', required=True, help="Sequence of the protein to be reverse-translated")
    parser.add_argument('--organism', '-o', default='Escherichia coli general', help="Organism for protein expression")
    args = parser.parse_args()

    # Call the function and print the result
    dna_seq = CodonTransformer_process(args.protein, args.organism)
    print(dna_seq)