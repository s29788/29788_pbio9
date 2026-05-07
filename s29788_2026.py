"""
s29788, 7.05.2026
program do generowania losowych sekwencji nukleotydowych z zapisem w formacie FASTA
"""

import random, math, string

nucleotides = ["A", "C", "G", "T"]

def generate_sequence(length: int) -> str:
    """
    Generuje losową sekwencję nukleotydową o podanej długości.
    Za pomocą metody random.randint() losuje pozycję nukleotydu w liście nucleotides.
    Nukleotyd o wylosowanej pozycji dokleja do sekwencji.

    Args:
        length (int): długość sekwencji nukleotydowej do wygenerowania
    Returns:
        str: losowo wygenerowana sekwencja nukleotydowa
    """
    seq = ""
    for i in (range(0, length)):
        rand = random.randint(0, len(nucleotides)-1)
        seq += nucleotides[rand]
    return seq

def calculate_stats(sequence: str) -> dict:
    """
    Oblicza procentową zawartość nukleotydów A, C, G, T oraz GC-content dla podanej sekwencji.
    Zwraca pusty słownik, jeśli sekwencja jest null.
    Usuwa imię z sekwencji za pomocą metody remove_name().
    Zlicza wystąpienia w sekwencji konkretnego nukleotydu za pomocą metody sequence.count(), otrzymaną liczbę przelicza na procent wystąpienia nukleotydu w sekwencji w porównaniu do liczby wszystkich nukleotydów w sekwencji.
    Oblicza GC-content przez dodanie wartości procentowych nukleotydów G i C.
    Formatuje za pomocą .2f przekazany parametr % (procentowe zawartości nukleotydów), żeby miał dwa miejsca po przecinku.

    Args:
        sequence (str): sekwencja DNA, dla której obliczane są statystyki

    Returns:
        dict: słownik zawierający statystyki procentowe
    """
    if len(sequence) == 0:
        return {}

    sequence = remove_name(sequence)

    a_count = sequence.count("A") / len(sequence) * 100
    c_count = sequence.count("C") / len(sequence) * 100
    g_count = sequence.count("G") / len(sequence) * 100
    t_count = sequence.count("T") / len(sequence) * 100
    gc_count = g_count + c_count
    a_stat = "%.2f" % a_count
    c_stat = "%.2f" % c_count
    g_stat = "%.2f" % g_count
    t_stat = "%.2f" % t_count
    gc_stat = "%.2f" % gc_count
    return {"A":str(a_stat)+"%", "C":str(c_stat)+"%", "G":str(g_stat)+"%", "T":str(t_stat)+"%", "GC-content":str(gc_stat)+"%"}

def insert_name(sequence: str, name: str) -> str:
    """
    Wstawia imię użytkownika w losowe miejsce sekwencji.
    Metoda lower() zamienia litery podane imię na małe, żeby można je było odróżnić od nukleotydów sekwencji.
    Zwraca samo imię, jeśli sekwencja jest null.
    Za pomocą metody random.randint() losuje pozycję nukleotydu w liście nucleotides.
    Dokleja imię po wylosowanej pozycji (sequence[:rand]) i dodaje resztę sekwencji

    Args:
        sequence (str): sekwencja, do której wstawiane jest imię
        name (str): imię wstawiane do sekwencji

    Returns:
        str: sekwencja z wstawionym imieniem
    """
    name.lower()

    if len(sequence) == 0:
        return name

    rand = random.randint(0, len(nucleotides))
    sequence = sequence[:rand] + name + sequence[rand:]
    return sequence

def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    """
    Formatuje sekwencję do standardu FASTA.
    Tworzy nagłówek odpowiedni dla formatu FASTA.
    Dokleja w linii pod nagłówkiem odpowiednio sformatowaną sekwencję DNA, tak, żeby w jednej linii sekwencji było tyle znaków, ile podano jako parametr.

    Args:
        seq_id (str): identyfikator sekwencji
        description (str): opis sekwencji
        sequence (str): sekwencja DNA
        line_width (int): maksymalna liczba znaków w linii

    Returns:
        str: sformatowany tekst FASTA
    """
    output = ">{0} {1}\n".format(seq_id, description)
    output += format_sequence(sequence, line_width)

    return output

def format_sequence(sequence: str, line_width: int) -> str:
    """
    Dzieli sekwencję na linie o podanej długości.
    Zwraca pusty string, jeśli podana sekwencja jest pusta.
    Oblicza liczbę linii, którą będzie mieć sformatowana sekwencja, dzieląc długość sekwencji przez liczbę znaków w jednej linii.
    Dla każdej linii dodaje odpowiednią ilość nukleotydów z podanej sekwencji, pobranych z fragmentu sekwencji,
    określonego przez pozycję start (początek sekwencji, dalej pierwszy nukleotyd po zakończeniu poprzedniej linii)
    i pozycję end (ilość znaków w jednej linii lub koniec sekwencji, jeśli to ostatnia linia, a jest krótsza od podanej liczby znaków w jednej linii).

    Args:
        sequence (str): sekwencja DNA
        line_width (int): liczba znaków w jednej linii

    Returns:
        str: sformatowana sekwencja
    """
    output = ""

    if len(sequence) == 0:
        return output

    start = 0
    if len(sequence) < line_width:
        end = len(sequence)
    else:
        end = line_width

    for i in range(0, math.ceil(len(sequence)/line_width)):
        for j in [sequence[start:end]]:
            output += j
        output += "\n"

        if start + line_width > len(sequence):
            start = len(sequence)
        else:
            start += line_width
        if end + line_width > len(sequence):
            end = len(sequence)
        else:
            end += line_width

    return output

def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    """
    Sprawdza, czy podana wartość jest liczbą całkowitą, jeśli wystąpi błąd podczas konwersji podanej wartości na int oraz czy podana wartość jest z określonego zakresu.

    Args:
        prompt (str): wartość podana przez użytkownika
        min_val (int): minimalna dozwolona wartość
        max_val (int): maksymalna dozwolona wartość

    Returns:
        int:
        1 - gdy dane są poprawne
        0 - gdy dane są niepoprawne
    """
    try:
        if int(prompt) >= 1 and int(prompt) <= 100000:
            return 1
        else:
            print("dlugosc sekwencji ma byc z zakresu [{0}, {1}]".format(min_val, max_val))
            return 0
    except ValueError:
        print("dlugosc sekwencji ma byc liczba calkowita")
        return 0

def check_for_white_spaces(input: str) -> int:
    """
    Sprawdza, czy napis zawiera białe znaki.

    Args:
        input (str): napis do sprawdzenia

    Returns:
        int:
            1 - brak białych znaków
            0 - napis zawiera spacje lub jest pusty
    """
    if input == "":
        return 0
    if " " in input:
        return 0
    else:
        return 1

def remove_name(sequence: str) -> str:
    """
    Usuwa z sekwencji małe litery odpowiadające imieniu.
    Zwraca pusty string, jeśli podana sekwencja jest pusta.
    Tworzy tabelę mapującą przy pomocy metody maketrans("", "", string.ascii_lowercase), która po użyciu metody translate() usuwa z sekwencji, na której jest wywołana małe litery.

    Args:
        sequence (str): sekwencja DNA

    Returns:
        str: sekwencja bez imienia
    """
    if len(sequence) == 0:
        return sequence

    table = str.maketrans("", "", string.ascii_lowercase)
    return sequence.translate(table)

def generate_compliment_sequence(sequence: str) -> str:
    """
    Generuje sekwencję komplementarną DNA.
    Tworzy tabelę mapującą przy pomocy metody maketrans({"A": "T", "T": "A", "C": "G", "G": "C"}), która po użyciu metody translate() zamienia w sekwencji, na której jest wywołana A<->T, G<->C.

    Args:
        sequence (str): oryginalna sekwencja DNA

    Returns:
        str: sekwencja komplementarna
    """
    sequence = remove_name(sequence)

    table = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C"})
    sequence = sequence.translate(table)

    return sequence

def generate_reverse_compliment_sequence(sequence: str) -> str:
    """
    Generuje odwrotnie komplementarną sekwencję DNA.
    Usuwa imię z podanej sekwencji za pomocą metody remove_name().
    Generuje komplementarną sekwencję DNA do podanej sekwencji za pomocą metody generate_compliment_sequence().
    Odwraca kolejność sekwencji przechodząc po niej wstecz (wartość -1 w [::-1])

    Args:
        sequence (str): oryginalna sekwencja DNA

    Returns:
        str: odwrotnie komplementarna sekwencja DNA
    """
    sequence = remove_name(sequence)
    sequence = generate_compliment_sequence(sequence)
    sequence = sequence[::-1]
    return sequence

def transcription(sequence: str) -> str:
    """
    Przeprowadza transkrypcję DNA.
    Zwraca pusty string, jeśli sekwencja jest pusta.
    Usuwa imię z sekwencji za pomocą metody remove_name().
    Zastępuje T przez U za pomocą metody replace().

    Args:
        sequence (str): sekwencja DNA

    Returns:
        str: sekwencja RNA
    """
    if len(sequence) == 0:
        return sequence

    seq = remove_name(sequence)
    return seq.replace("T", "U")

def main():
    """
    Główna funkcja programu odpowiedzialna za:
    - pobieranie danych od użytkownika,
    - generowanie sekwencji,
    - zapis do pliku FASTA,
    - wyświetlanie statystyk,
    - generowanie i zapisywanie do pliku sekwencji komplementarnej i odwrotnie komplementarnej,
    - transkrypcję sekwencji i zapisywanie wyniku do pliku
    """

    var = 0
    seq_length = 0
    while var != 1:
        seq_length = str(input("podaj dlugosc sekwencji "))
        var = validate_positive_int(seq_length)
    var = 0
    while var != 1:
        seq_id = str(input("podaj id sekwencji "))
        var = check_for_white_spaces(seq_id)
        if var == 0:
            print("id nie moze zawierac bialych znakow")

    seq_desc = str(input("podaj opis sekwencji "))

    seq = generate_sequence(int(seq_length))

    name = str(input("podaj imie "))
    seq = insert_name(seq, name)

    file_name = "{0}.fasta".format(seq_id)
    with open(file_name, mode="w", encoding="utf-8") as file:
        file.write(format_fasta(seq_id, seq_desc, seq)+"\n")
    print("zapisano sekwencje do pliku " + file_name)

    print(calculate_stats(seq))

    with open(file_name, mode="a", encoding="utf-8") as file:
        file.write("sekwencja komplementarna\n" + format_sequence(generate_compliment_sequence(seq), 80)+"\n")
        file.write("sekwencja odwrotnie komplementarna\n" + format_sequence(generate_reverse_compliment_sequence(seq), 80)+"\n")

    with open(file_name, mode="a", encoding="utf-8") as file:
        file.write("po transkrypcji\n" + format_sequence(transcription(seq), 80) + "\n")

if __name__ == "__main__":
    main()
