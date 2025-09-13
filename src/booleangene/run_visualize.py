import os
from .modules.visualization_helper import visualize

def main(input: str, output: str, labels: bool=False, seed=None) -> None:
    
    output_path = os.path.expanduser(output)
    fig_orig, fig_scrambled = visualize(input, labels, seed)
    fig_orig.savefig(os.path.join(output_path, "network_ORIGINAL.png"), dpi=400, bbox_inches='tight')
    fig_scrambled.savefig(os.path.join(output_path, "network_SCRAMBLED.png"), dpi=400, bbox_inches='tight')
    return