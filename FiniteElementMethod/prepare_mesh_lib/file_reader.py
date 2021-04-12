import re
import os


class FileReader:
    def __init__(self, filename, mesh_dir, raw_mesh_dir,
                 new_dir='new_mesh', dim=2):
        # Имя файла который был создан фидесисом,
        # где хранится сетка.
        self.filename = filename
        self.mesh_dir = mesh_dir
        self.raw_mesh_dir = raw_mesh_dir
        self.new_dir = new_dir
        self.dim = dim

    def parse_nodes(self):
        # Требуется распарсить и записать Nodes
        # в csv файл, удобный для чтения в c++.

        filename = self.raw_mesh_dir + '/' + self.filename
        nodes = []

        print('Parsing file for nodes start...')

        with open(filename, 'r') as mesh:
            for line in mesh:

                # Ищем строчку с которой начинаются Nodes,
                # сохраняем все следующие строчки, пока
                # не встретим '$' - символ начала следующего блока.
                if line.strip() == '*NODE':
                    line = mesh.readline()[1:]  # Т.к. там вначале стоит символ $.
                    line = mesh.readline()  # Чтобы сразу считывать числа.
                    while line.strip() != '$':
                        nodes.append(line.strip())
                        line = mesh.readline()

        print('Parsing file for nodes end.')

        self.raw_nodes = nodes
        print('OK')

    def prepare_nodes(self):
        # Требуется преобразовать массив raw_nodes, состоящий
        # из строчек, в каждой из которой координаты точек

        try:
            nodes = self.raw_nodes
        except Exception as e:
            print('Use method parse_nodes first!!')
            return

        print('Prepare and writing nodes in new file start...')

        with open(self.dir_name + '/nodes.txt', 'w') as write_file:
            #print('Attention dimension of mesh == {0}'.format(self.dim))
            for node in nodes:
                line = ' '.join(re.split('\s+', node)[:self.dim+1])
                print(line)
                print(type(line))
                print(len(line))
                write_file.write(line + '\n')

        print('Prepare and writing nodes in new file end.')
        print('OK')

    def parse_elements(self):
        # Требуется распарсить и записать Elements
        # в csv файл, удобный для чтения в c++.

        filename = self.raw_mesh_dir + '/' + self.filename
        elements = []

        print('Parsing file for elements start...')

        with open(filename, 'r') as mesh:
            for line in mesh:

                # Ищем строчку с которой начинаются Elements,
                # сохраняем все следующие строчки, пока
                # не встретим '*END' - символ начала следующего блока.
                if line.strip() == '*ELEMENT_SHELL':

                    line = mesh.readline()
                    while line.strip() != '*END':
                        elements.append(line.strip())
                        line = mesh.readline()

        print('Parsing file for elements end.')

        self.raw_elements = elements
        print('OK')

    def prepare_elements(self):
        # Требуется преобразовать массив raw_elements, состоящий
        # из строчек, в каждой из которой номера узлов,
        # в удобный для c++ csv файл, где эти данные через запятую.

        try:
            elements = self.raw_elements
        except Exception as e:
            print('Use mehtod parse_elements first!!')
            return

        print('Prepare and writing elements in new file start...')

        with open(self.dir_name + '/elements.txt', 'w') as write_file:
            #print('Attention dimension of mesh == {0}'.format(self.dim))
            for element in elements:
                line = ' '.join(re.split('\s+', element))
                write_file.write(line + '\n')

        print('Prepare and writing elements in new file end.')
        print('OK')

    def make_directory(self):
        # Делаю новую директорию, куда
        # буду сбрасывать данные связанные с этой сеткой.
        # Требуется начинать работу с этой команды.
        print('Make new directory named: ' + self.new_dir
              + ' in ' + self.mesh_dir)

        self.dir_name = self.mesh_dir + '/' + self.new_dir
        try:
            os.mkdir(self.dir_name)
        except Exception as e:
            print('This directory already exist!')
            return
        print('OK')

    def make_good(self):
        # Переносит исходный файл в новую директорию.
        # TODO: добавить описание сетки, типа README.
        # Требуется заканчивать работу этой командой.

        print('Move base file in new directory.')
        try:
            os.replace(self.raw_mesh_dir + '/' + self.filename,
                       self.dir_name + '/' + self.filename)
        except Exception as e:
            print('Some troubles with moving ' + self.filename
                  + ' to ' + self.new_dir)
