import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter molecules with PAINS")
    subparsers = parser.add_subparsers(dest= 'subcmd', help='subcommands', metavar='SUBCOMMAND')
    subparsers.required = True

    


