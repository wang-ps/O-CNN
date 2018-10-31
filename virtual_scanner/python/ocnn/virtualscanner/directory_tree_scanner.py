try:
    from Queue import Queue
except ModuleNotFoundError:
    from queue import Queue

import os

from ocnn.virtualscanner import VirtualScanner
from threading import Thread

class DirectoryTreeScanner:
    """
        Walks a directory and converts off/obj files to points files.
    """
    def __init__(self, view_num=6, flags=False, normalize=False):
        """
            Initializes DirectoryTreeScanner

            Args:
                view_num (int): The number of view points to scan from.
                flags (bool): Indicate whether to ouput normal flipping flag.
                normalize (bool): Normalize maximum extents of mesh to 1.
        """
        self.view_num = view_num
        self.flags = flags
        self.normalize = normalize
        self.scan_queue = Queue()

    def _scan(self):
        """
            Creates VirtualScanner object and creates points file from obj/off
        """
        while(True):
            input_path, output_path = self.scan_queue.get()

            print('Scanning {0}'.format(input_path))
            scanner = VirtualScanner(input_path,
                                     self.view_num,
                                     self.flags,
                                     self.normalize)
            scanner.save(output_path)
            self.scan_queue.task_done()

    def scan_tree(self, input_base_folder, output_base_folder='', num_threads=1):
        """
            Walks directory looking for obj/off files. Outputs points files for
            found obj/off files.

            Args:
                input_base_folder (str): Base folder to scan
                output_base_folder (str): Base folder to output points files in
                    mirrored directory structure. If not specified will output
                    to input_base folder.
                num_threads (int): Number of threads to use to convert obj/off
                    to points
        """
        if not output_base_folder:
            output_base_folder = input_base_folder

        if not os.path.exists(output_base_folder):
            os.mkdir(output_base_folder)

        for i in range(num_threads):
            scan_thread = Thread(target=self._scan)
            scan_thread.daemon = True
            scan_thread.start()

        for root, _, files in os.walk(input_base_folder):
            rel_path = os.path.relpath(root, input_base_folder)

            output_folder = os.path.join(output_base_folder, rel_path)
            if not os.path.exists(output_folder):
                os.mkdir(output_folder)

            for filename in files:
                basename, extension = os.path.splitext(filename)
                extension = extension.lower()
                if extension == '.obj' or extension == '.off':
                    outfilename = basename + '.points'
                    input_path = os.path.join(root, filename)
                    output_path = os.path.join(output_folder, outfilename)
                    self.scan_queue.put((input_path, output_path))
        self.scan_queue.join()
