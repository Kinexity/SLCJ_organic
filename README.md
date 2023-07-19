# SLCJ_organic

Instalacja tak jak każdego projektu cmake.

Wymagania:
- C++20
- boost z modułem program_options
- GEANT4 10.5 (taką wersję testowałem)

# Uruchamianie

Lista opcji wiersza polecenia z przykładowymi wartościami:

--skipIfDataExists - flaga, która powoduje, że program pominie symulację jeżeli wykonał już jedną z takimi parametrami

--dataOverwrite - flaga, która powoduje, że dane z symulacji z indeksem 0 są nadpisywane

--isGeantino - flaga wywołująca symulację z cząstką geantino. Zapisuje geometrię detektora w pliku geantino.txt w formacie X Y Z R G B. Polecam program MeshLab do otwierania

--foodVolume 10 - opcja ustalająca objętość pożywki (w cm3)

--cutValue 0.01 - opcja decydująca o parametrze cut value symulacji, który decyduje o dokładności (w mm). Wyższa dokładność = wolniejsza symulacja. Polecam nie używać ponieważ domyślna wartość powinna być ok. Typowy zakres: 0.01-1 mm

--energy 9 - energia wiązki (w MeV). Pomijane dla geantino

--physicsListName empenelope - physics list decydujący o zakresie symulacji. Dostępne: local, emlivermore, empenelope. Polecam pomijać.

--numberOfEvent 1000000 - liczba cząstek symulowanych sekwencyjnie w przebiegu symulacji.

Wszystkie opcje są opcjonalne. Jeżeli jakaś nie jest przekazana do programu to użyje on wartości domyślnej.

Przykładowe wywołanie programu:

.\SLCJ.exe --numberOfEvent 1000000 --cutValue 0.01 --foodVolume 10

# Dane wyjściowe symulacji

Program tworzy folder results_SLCJ pod ścieżką, z której jest uruchamiany i tam tworzy odpowiedni folder dla danych z symulacji np. event_G4EmLivermore_0.01mm_10cm3_9MeV_5. W takim przykładowym folderze będą następujące pliki:
- elements.csv - lista zaimplementowanych pierwiastków
- isotopes.csv - lista zaimplementowanych izotopów
- materials.csv - lista zaimplementowanych materiałów
- volumes.csv - lista zaimplementowanych objętości w symulowanej strukturze
- metadata.bin - lista początkowych parametrów wysyłanych cząstek (nieistotne więc pomijam wytłumaczenie formatu tych danych)
- event_G4EmLivermore_0.01mm_10cm3_9MeV_energyDeposit.bin - plik binarny zawierający wszystkie poszczególne punkty depozycji energii w objętości napromieniowywanych komórek. W pliku wszystkie wartości to liczby są typu 64-bit float odpowiednio w kolejności (energia,x,y,z), więc w pliku mamy po kolei wartości E1,x1,y1,z1,E2,x2,y2,z2,E3,x3,y3,z3 i tak dalej.

Tutaj jest implementacja w pythonie prostego programu, który otwiera te pliki z danymi i tworzy proste wizualizacje (wymaga zmiany ścieżek i skonfigurowania pod własne dane): https://github.com/Kinexity/SLCJ_dataAnalysis/blob/master/SLCJ_dataAnalysis/SLCJ_dataAnalysis.py
