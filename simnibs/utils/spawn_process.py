# import sys
# from threading import Thread
# from multiprocessing import Process
# from queue import Queue
import logging
from subprocess import Popen, PIPE, STDOUT, CalledProcessError

from .simnibs_logger import logger

# def spawn_process(cmd, new_thread=False, lvl=logging.INFO, shell=False,
#                   new_process=False, check=True):
#     """Spawn a new process and communicate its output.
    
#     PARAMETERS
#     ----------
#     cmd : str
#         The command to be executed.
#     new_thread : bool, optional
#         By default, the child process blocks the main process while running,
#         however, by setting this option to true the child progress will be
#         spawned in a new thread, hence not blocking the main process (default =
#         False).
#     lvl : int, optional
#         Logging level. Only applies if  new_thread == False (default = logging.INFO).
#     shell: bool, optional
#         Whether to use shell mode. Default: False (forced to be True on Windows)
#     new_process: bool, optional
#         Starts command in a new process. Default: False
#     ckeck: bool, optional
#         If new_thread=False and new_thread=False, checks the output and raises an error
#         if the reurncode is not 0. Default: True
#     RETURNS
#     ----------
#     exit_status : int
#         Return code of the process.
#     """
#     logger.debug(f"Running {cmd}")
#     if new_thread or new_process:
#         ON_POSIX = "posix" in sys.builtin_module_names

#         def enqueue_output(out, queue):
#             for line in iter(out.readline, b''):
#                 try:
#                     queue.put(line.decode())
#                 except UnicodeDecodeError:
#                     queue.put('Could not print line')
#             out.close()

#         p = Popen(cmd, stdout=PIPE,
#                   bufsize=1, close_fds=ON_POSIX,
#                   shell=True)

#         q = Queue()
#         if new_thread:
#             t = Thread(target=enqueue_output, args=(p.stdout, q))
#         if new_process:
#             t = Process(target=enqueue_output, args=(p.stdout, q))
#         t.daemon = True  # thread dies with the program
#         t.start()
#         return t

#     else:
#         p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
#         for line in iter(p.stdout.readline, b''):
#             # to prevent massive output to the logfile retain only the last
#             # line starting with \r
#             try:
#                 line = line.decode().split("\r")[-1]
#             except UnicodeDecodeError:
#                 logger.debug('Could not print line')
#                 continue

#             # remove one \n since logger.log adds one itself
#             if line[-1] == "\n":
#                 line = line[:-1]

#             logger.log(lvl, line)

#         p.wait()
#         if check and p.returncode != 0:
#             raise CalledProcessError(
#                 p.returncode, cmd
#             )

#         return p.returncode


def spawn_process(cmd, lvl=logging.INFO):
    """Spawn a new process and log its output.
    
    PARAMETERS
    ----------
    cmd : list
        The command to be executed, passed as a list
    lvl : int, optional
        Logging level (default = logging.INFO).

    """
    logger.debug(f"Running {cmd}")
    p = Popen(cmd, stdout=PIPE, stderr=STDOUT)
    for line in iter(p.stdout.readline, b''):        
        # get all output
        line = line.decode('ASCII', errors='ignore').replace('\r', '')
        # remove one \n since logger.log adds one itself
        if line[-1]=="\n":
             line = line[:-1]
        logger.log(lvl, line)
    p.wait()
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, cmd)
    return 0
