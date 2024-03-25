import yaml
import os

def test_config():
    # Check if the config file exists
    assert os.path.exists('config.yaml'), "config.yaml does not exist"

    # Load the config file
    with open('config.yaml') as f:
        config = yaml.safe_load(f)

    # Check for essential keys in the config file
    assert 'samples' in config, "Config file missing 'samples' key"
    assert 'reference' in config, "Config file missing 'reference' key"

    # You can also check if the reference file exists
    assert os.path.exists(config['reference']), f"Reference file {config['reference']} does not exist"

    # If you want, check for the existence of sample files
    for sample, paths in config['samples'].items():
        for path in paths.values():
            assert os.path.exists(path), f"Sample file {path} for sample {sample} does not exist"
