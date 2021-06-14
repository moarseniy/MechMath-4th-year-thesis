# Install script for directory: C:/Users/mokin/Desktop/cuda-template-master3D/cuda-template-master

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/FemProject")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("C:/Users/mokin/Desktop/cuda-template-master3D/build-cuda-template-master-Desktop_x86_windows_msvc2019_pe_64bit-Release/Linal2/cmake_install.cmake")
  include("C:/Users/mokin/Desktop/cuda-template-master3D/build-cuda-template-master-Desktop_x86_windows_msvc2019_pe_64bit-Release/CudaLib/cmake_install.cmake")
  include("C:/Users/mokin/Desktop/cuda-template-master3D/build-cuda-template-master-Desktop_x86_windows_msvc2019_pe_64bit-Release/fem/cmake_install.cmake")
  include("C:/Users/mokin/Desktop/cuda-template-master3D/build-cuda-template-master-Desktop_x86_windows_msvc2019_pe_64bit-Release/VTKfile/cmake_install.cmake")
  include("C:/Users/mokin/Desktop/cuda-template-master3D/build-cuda-template-master-Desktop_x86_windows_msvc2019_pe_64bit-Release/femfunc/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "C:/Users/mokin/Desktop/cuda-template-master3D/build-cuda-template-master-Desktop_x86_windows_msvc2019_pe_64bit-Release/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
