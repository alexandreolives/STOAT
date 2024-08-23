import time

def determine_str(s: str, length_s : int, i: int) -> tuple[int, int]:
    """Extract an integer from a string starting at index i."""
    start_idx = i
    while i < length_s and s[i] not in ['>', '<']:
        i += 1
    return i, s[start_idx:i]

def decompose_string(s: str) :
    """Decompose a string with snarl information."""
    print("s : ", s)
    result = []
    i = 0
    length_s = len(s)
    prev_int = None
    prev_sym = None
    
    while i < length_s:
        start_sym = s[i]
        i += 1
        i, current_int = determine_str(s, length_s, i)

        if prev_int is not None and prev_sym is not None:
            result.append(f"{prev_sym}{prev_int}{start_sym}{current_int}")
        
        prev_int = current_int
        prev_sym = start_sym
    
    print(result)
    return result

string = ">1322<1323>1323<1323>1325"

start = time.time()
snarl_decomposed = decompose_string(string)
print(f"Time spend : {time.time() - start}")

import pandas as pd

data = {
    '>1>2>7': [0, 0, 0, 0],
    '>1>2>3>5>6': [1, 1, 0, 0]
}
index = ['SRR933569', 'SRR933570', 'SRR933571', 'SRR933572']
df = pd.DataFrame(data, index=index)
print(df)

import pandas as pd

# 1. Create the initial DataFrame with boolean values
data = {
    '>1>2>7': [0, 0, 0, 0],
    '>1>2>3>5>6': [1, 1, 0, 0]
}
index = ['SRR018517', 'SRR018521', 'SRR018574', 'SRR018580']
df = pd.DataFrame(data, index=index)

# Display the initial DataFrame
print("Initial DataFrame:")
print(df)

# 2. Pheno dictionary with target values
pheno = {
    'SRR933569': 47.52, 
    'SRR933570': 75.1020408163265, 
    'SRR933571': 63.68, 
    'SRR933572': 66.72, 
    'SRR933573': 42.3829787234043, 
    'SRR933575': 55.52, 
    'SRR933577': 53.2525252525253, 
    'SRR933578': 68.8, 
    'SRR933579': 44.32, 
    'SRR933580': 50.88
}

# 3. Map the target values from the dictionary to the DataFrame based on index
df['Target'] = df.index.map(pheno)

# Display the updated DataFrame
print("\nUpdated DataFrame with Target column:")
print(df)