import subprocess
import shlex
import sys
from threading import Thread
try:
    from Queue import Queue, Empty
except:
    from queue import Queue, Empty


from .simnibs_logger import logger


def run_command(command, realtime_output=False):
    """ Run a command and logs it

    Parameters
    --------------------------------------
    command: list or string
        list of strings with arguments, subprocess style
        or a string to be split and transformed into a list of arguments

    realtime_output (optional): bool
        True if output is to be returned in realtime. Default = False

    Returns:
    --------------------------------------
    str
        Process STDOUT and STDERR output

    Raises
    ----------------------------------------
    OSError
        if there was a problem running the comand
    """
    if isinstance(command, str):
        args = shlex.split(command)
    else:
        args = command

    logger.info('Executing: \n' + ' '.join(args))

    try:
        command_line_process = subprocess.Popen(
            ' '.join(args),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)

        if realtime_output:
            process_output = ''
            while True:
                rc = command_line_process.poll()
                output = command_line_process.stdout.readline().decode('ascii')
                if output != '':
                    logger.info(output.rstrip('\n'))
                    process_output += output

                if output == '' and rc is not None:
                    break

        else:
            process_output, _ = command_line_process.communicate()
            rc = command_line_process.returncode
            logger.info(process_output)

    except OSError:
        logger.error('Could not execute command')
        raise OSError('Could not execute command:\n' + ' '.join(args))

    if rc == 0:
        logger.info('Execution finished')

    else:
        logger.error('Error while executing command:\n' + ' '.join(args))
        logger.error('Command Output:\n' + process_output)
        raise OSError('Error executing command:\n' + ' '.join(args))

    return process_output

def run_command_new_thread(command, daemon=False):
    if isinstance(command, str):
        args = shlex.split(command)
    else:
        args = command

    logger.info('Executing: \n' + ' '.join(args))
    ON_POSIX = "posix" in sys.builtin_module_names

    def enqueue_output(out, queue):
        for line in iter(out.readline, b''):
            queue.put(line)
        out.close()

    p = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        bufsize=1, close_fds=ON_POSIX)
    q = Queue()
    t = Thread(target=enqueue_output, args=(p.stdout, q))
    t.daemon = daemon  # thread dies with the program
    t.start()
