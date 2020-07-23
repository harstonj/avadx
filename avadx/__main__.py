#!/usr/bin/env python3

import sys
if sys.version_info[0] == 3 and sys.version_info[1] >= 8:
    pass
else:
    raise Exception("mi-faser requires Python >= 3.8. This version is: %s" % ".".join('%s' % _ for _ in sys.version_info[:3]))


def main():
    from .pipeline import init
    init()


if __name__ == "__main__":
    main()
