"""Szacowanie liczby pi metodą Monte Carlo."""

import random


def oszacuj_pi(liczba_losowan: int) -> float:
    punkty_w_kole = 0

    for _ in range(liczba_losowan):
        x = random.random()
        y = random.random()

        if x * x + y * y <= 1:
            punkty_w_kole += 1

    return 4 * punkty_w_kole / liczba_losowan


liczba_losowan = 1_000_000
wynik = oszacuj_pi(liczba_losowan)

print(f"Oszacowanie liczby pi: {wynik}")