import ijson
import sys


def main():
    input_json_file = sys.argv[1]

    with open(input_json_file, "r") as input_json:
        ijson_generator = ijson.items(input_json, "genes.item")
        for gene in ijson_generator:
            for key in gene:
                print(gene["id"], key, type(gene[key]), sep="\t")


if __name__ == "__main__":
    main()
