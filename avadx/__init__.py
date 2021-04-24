from pathlib import Path


name = "avadx"
__author__ = 'mmiller'
__version__ = '2.0.11'
__releasedate__ = '04/24/21'
__build__ = {
    line.strip().split()[-1]:
    line.strip().split()[0] for line in ((Path(__file__).parent / '.build').open() if (Path(__file__).parent / '.build').exists() else [])
}.get('HEAD', None)
__all__ = [
    'pipeline',
    'helper',
    'logger',
]
