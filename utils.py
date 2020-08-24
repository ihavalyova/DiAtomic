import os
import shutil
import datetime


class Utils:
    def __init__(self):
        pass

    @classmethod
    def createBackup(cls, ref_file):

        backup_folder = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), '.backup'
        )

        os.makedirs(backup_folder, exist_ok=True)

        fname, fext = os.path.splitext(ref_file)

        name_bkp = '_'.join([fname, 'bkp', cls.getDatetime()]) + fext

        shutil.copy2(ref_file, os.path.join(backup_folder, name_bkp))

    @classmethod
    def getDatetime(cls):

        now = datetime.datetime.now()
        # milli = str(int(round(now.microsecond / 1.0e3)))  # us to ms

        timelist = [
            now.year, now.month, now.day, now.hour,
            now.minute, now.second, now.microsecond
        ]

        return '_'.join(str(t) for t in timelist)

    @classmethod
    def get_current_dir(cls):

        # will not work if the program is used as a module
        # return os.path.abspath(os.path.dirname(__file__))

        return os.getcwd()

    # @classmethod
    # def calculate_rms(cls, yobs, ycal):
    #     s = np.sum(np.square(yobs-ycal)) / yobs.shape[0]
    #     return math.sqrt(s)
    # @classmethod
    # def calculate_dimless_rms(cls, yobs, ycal, yunc, weighting):
    #     diff_square = np.square(yobs - ycal)
    #     if not weighting:
    #         weights = 1.0 / np.square(yunc)
    #     else:
    #         weights = 1.0 / (np.square(yunc) + (diff_square / 3.0))

    #     s = np.sum(diff_square * weights) / yobs.shape[0]

    #     return math.sqrt(s)
