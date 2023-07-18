# SLCJ_organic

Instalacja tak jak każdego projektu cmake.

Wymagania:
- boost z modułem program_options
- GEANT4 10.5 (taką wersję testowałem)

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
