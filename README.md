# Описание
* Файл `rational.h` содержит реализацию рациональных чисел(они почти нигде не используются)
* Файл `polynom.h` содержит реализацию многочлена(Для подсчета хар. многочлена)
* Файл `matrix.h` содержит реализацию матрицы и ее функций
* В файле `check.cpp` содержится `namespace` с тестированием(он не нужен для дз), а также с нужными вычислениями(по каждой задаче поотдельности)
# Как компилировал я
* Я компилировал `check.cpp` только под `C++17`: 
```
g++-12 ./check.cpp -std=c++17 -Wall -o ./check
```
* Warning'и есть, но их всего 3, и не очень серьезные, так что можно забить.