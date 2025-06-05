""" Select new FASTQ files from gcloud bucket """

import os

with open('fastq_urls.txt', 'r') as f:
    fastq_urls = f.read().splitlines()

fastq_urls = [url for url in fastq_urls if '.fastq' in url]

sample_ids = list(set([url.split('/')[-1].split('_R')[0] for url in fastq_urls]))
print(f'Found {len(sample_ids)} samples in gcloud bucket')

completed_samples = [file.split('.covariants.tsv')[0] for file in os.listdir('covariants')]
print(f'Found {len(completed_samples)} completed samples')

new_samples = [sample for sample in sample_ids if sample not in completed_samples]
print(f'Found {len(new_samples)} new samples')

new_fastq_urls = [url for url in fastq_urls if url.split('/')[-1].split('_R')[0] in new_samples]

with open('new_fastq_urls.txt', 'w') as f:
    for url in new_fastq_urls:
        f.write(f'{url}\n')