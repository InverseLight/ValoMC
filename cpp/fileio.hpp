#ifndef _FILEIO_H_
#define _FILEIO_H_

#include <sstream>

template <typename F>
void readAndResize(FILE *file, int height, int width, bool halt, std::vector<Array<F> *> items, std::vector<std::string> names)
{
  bool success = true;
  char line[5012];

  // read the name tag and rewind if it does not match
  fpos_t position;
  fgetpos (file, &position);
  
  bool should_rewind = false;
  
  if(fgets(line, 5011, file) == 0) {
       return;
  };

  std::stringstream linestream;
  linestream << line;
  std::istringstream iss();
  std::vector<std::string> infoline;
  std::string word;

  while (linestream >> word)
  {
    infoline.push_back(word);
  }

  if(infoline.size() == names.size() ) {
     for (unsigned int i = 0; i < items.size(); i++) {
      	if(infoline[i].compare(names[i])) {
           should_rewind = true;
      	}
     }
  } else {
     should_rewind = true;
  }
  
  if(should_rewind) {
     fsetpos(file, &position);
     return;
  }

  for (unsigned int i = 0; i < items.size(); i++)
  {
    items[i]->resize(height, width);
  }
  for (int_fast64_t i = 0; i < height; i++)
  {
    if (!feof(file) && fgets(line, 5011, file) != NULL)
    {
      std::stringstream linestream;
      linestream << line;
      try
      {
        for (int_fast64_t j = 0; j < items.size(); j++)
        {
          for (int_fast64_t k = 0; k < width; k++)
          {
            F *entry = &(*items[(int)j])(i, k);
            linestream >> *entry;
            //  std::cout << names[j] << "[" << i << "," << k << "]"
            //            << ":" << std::setprecision(30) << *entry << ".";
          }
          //  std::cout << "\n";
        }
      }
      catch (...)
      {
        success = false;
        break;
      }
    }
    else
    {
      success = false;
      break;
    }
  }
  if (!success)
  {
    for (unsigned int i = 0; i < items.size(); i++)
    {
      items[i]->destroy();
    }
  }
  else
  {
    for (unsigned int i = 0; i < items.size(); i++)
    {
      std::cout << std::right << std::setw(12) << names[i].c_str();
      std::cout << "   "
                << "(" << items[i]->Nx << "x" << items[i]->Ny << ")\n";
    }
  }
  if (halt && !success)
  {
    printf("Reading error!\n");
    // todo throw error
    exit(1);
  }
}

template <typename F>
void readAndResize(FILE *file, int height, int width, bool halt, Array<F> *array, std::string name)
{
  std::vector<Array<F> *> items;
  std::vector<std::string> names;
  items.push_back(array);
  names.push_back(name);
  readAndResize(file, height, width, halt, items, names);
}

template <typename F>
void readAndResize(FILE *file, int height, int width, bool halt, Array<F> *array1, Array<F> *array2, Array<F> *array3, Array<F> *array4, std::string name1, std::string name2, std::string name3, std::string name4)
{
  std::vector<Array<F> *> items;
  std::vector<std::string> names;
  items.push_back(array1);
  names.push_back(name1);
  items.push_back(array2);
  names.push_back(name2);
  items.push_back(array3);
  names.push_back(name3);
  items.push_back(array4);
  names.push_back(name4);
  readAndResize(file, height, width, halt, items, names);
}

#endif
