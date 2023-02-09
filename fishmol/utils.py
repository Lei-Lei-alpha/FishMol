from IPython.display import clear_output


def update_progress(progress):
  
    bar_length = 20
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1
        
    block = int(round(bar_length * progress))

    clear_output(wait = True)
    text = "Progress: [{0}] {1:.1f}%".format( "■" * block + "○" * (bar_length - block), progress * 100)
    print(text)

# def to_ase_atoms()

def to_sublists(lst, length=2):
    """
    Split a list to sublists with specified length: e.g. a = [a,b,c,d,e]
    to_sublists(a) => [[a,b], [b,c], [c,d], [d,e]] 
    """
    return [lst[i:i+length] for i in range(len(lst)+1-length)]
