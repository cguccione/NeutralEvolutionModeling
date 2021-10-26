import click
from neuevo import __version__
from importlib import import_module


def _terribly_handle_brokenpipeerror():
    # based off http://stackoverflow.com/a/34299346
    import os
    import sys
    sys.stdout = os.fdopen(1, 'w')


@click.group()
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx):
    ctx.call_on_close(_terribly_handle_brokenpipeerror)


import_module('neuevo.scripts._neufit')