import os
import subprocess


def create_recovery_file(fetch_directory):
    """Copia fetch.txt a recovery.txt y borra fetch.txt. Devuelve la ruta de recovery.txt."""

    fetch_file = os.path.join(fetch_directory, "fetch.txt")
    recovery_file = os.path.join(fetch_directory, "recovery.txt")

    if not os.path.isfile(fetch_file):
        print("no fetch file")
        return None

    with open(fetch_file, "r") as f_fetch, open(recovery_file, "w") as f_recovery:
        for line in f_fetch:
            f_recovery.write(line)

    os.remove(fetch_file)
    return recovery_file


def splitting_fetch_file(recovery_file, n_split=100):
    """Generador: devuelve listas de n_split líneas de recovery.txt."""

    chunk = []
    with open(recovery_file, "r") as f:
        for line in f:
            chunk.append(line)
            if len(chunk) == n_split:
                yield chunk
                chunk = []
    if chunk:  # último bloque, más corto
        yield chunk


def processing(directory, n_split=100, max_workers=10):
    fetch_directory = os.path.join(directory, "ncbi_dataset")
    recovery_file = create_recovery_file(fetch_directory)
    if recovery_file is None:
        return

    new_fetch_file = os.path.join(fetch_directory, "fetch.txt")

    for i, chunk in enumerate(splitting_fetch_file(recovery_file, n_split)):
        # fetch.txt solo con este bloque
        with open(new_fetch_file, "w") as f:
            f.writelines(chunk)

        print(f"Bloque {i}: {len(chunk)} genomas")
        subprocess.run(
            ["datasets", "rehydrate", "--directory", directory,
             "--max-workers", str(max_workers)],
            check=True,
        )

    # al terminar, restaurar el fetch.txt completo
    os.replace(recovery_file, new_fetch_file)