from os import environ
import logging
import inspect


LOG_LEVEL = logging.INFO
LOG_FILE = None


class Logger:

    def __init__(self, name, level=None):
        level_checked = environ.get("LOGLEVEL", LOG_LEVEL) if level is None else level
        self.level = level_checked if type(level_checked) == str else getLevelName(level_checked)
        self.log = logging.getLogger(name)
        self.log.setLevel(self.level)

    def addConsoleHandler(self, level=None, formatted=True):
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level if level else self.level)
        if formatted and console_handler.level > 0:
            if console_handler.level >= 20:
                formatter = logging.Formatter(
                    "[ %(levelname)8s ] --- %(message)s"
                )
            else:
                formatter = logging.Formatter(
                    "%(asctime)s - %(filename)15s:%(lineno)-4s - [ %(levelname)8s ] - %(funcName)15s --- %(message)s"
                )
            console_handler.setFormatter(formatter)
        self.log.addHandler(console_handler)

    def addFileHandler(self, filename, level=None):
        file_handler = logging.FileHandler(filename)
        file_handler.setLevel(level if level else self.level)
        if file_handler.level >= 20:
            formatter = logging.Formatter(
                "[ %(levelname)8s ] --- %(message)s"
            )
        else:
            formatter = logging.Formatter(
                "%(asctime)s - %(filename)15s:%(lineno)-4s - [ %(levelname)8s ] - %(funcName)15s --- %(message)s"
            )
        file_handler.setFormatter(formatter)
        self.log.addHandler(file_handler)

    def addLoggers(self):
        if LOG_FILE:
            self.addConsoleHandler(logging.ERROR)
            self.addFileHandler(LOG_FILE)
        else:
            self.addConsoleHandler()

    def getFrame(self):
        return inspect.currentframe().f_back.f_code

    def getLogger(self):
        return self.log

    def getRootLogger(self):
        return logging.getLogger()


def update_logger(log_level, log_file):
    global LOG_LEVEL, LOG_FILE
    LOG_LEVEL = log_level
    LOG_FILE = log_file


def getLevelName(level):
    return logging.getLevelName(level)
