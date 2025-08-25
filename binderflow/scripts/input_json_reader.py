#! /usr/bin/env python3

import json
import argparse

def read_json_file(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    for key, value in data.items():
        # Ensure values with spaces are enclosed in quotes
        print(f'{key}="{value}"')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read and print JSON file contents.")
    parser.add_argument("file_path", type=str, help="Path to the JSON file to read.")
    args = parser.parse_args()

    file_path = args.file_path
    read_json_file(file_path)