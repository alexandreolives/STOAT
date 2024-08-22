import time
import re

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
