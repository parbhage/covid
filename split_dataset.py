import pandas as pd
from sklearn.model_selection import train_test_split

# 1. Load dataset
file_path = "/home/salix/Documents/covid_new/covid_R.txt"   #change to your dataset path
df = pd.read_csv(file_path, sep=",")

# 2. Train-test split (stratified by 'class')
train_df, test_df = train_test_split(
    df,
    test_size=0.2,
    random_state=42,
    stratify=df['Class']
)

# 3. Save to your desired folder
output_dir = "/home/salix/Documents/covid_new/" # change to your output folder
train_csv_path = output_dir + "train_dataset.csv"   # change to your desired train file name
test_csv_path = output_dir + "test_dataset.csv" # change to your desired test file name
train_txt_path = output_dir + "train_dataset.txt"       
test_txt_path = output_dir + "test_dataset.txt"

train_df.to_csv(train_csv_path, index=False)
test_df.to_csv(test_csv_path, index=False)
train_df.to_csv(train_txt_path, sep="\t", index=False)
test_df.to_csv(test_txt_path, sep="\t", index=False)

print("✅ Train file saved to:", train_csv_path)
print("✅ Test file saved to:", test_csv_path)
