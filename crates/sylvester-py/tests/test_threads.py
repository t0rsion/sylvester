"""The GIL is released around every computation."""

import threading

import sylvester


def test_python_runs_while_a_basis_computes(cyclic):
    """The worker computes cyclic-8 under a one second deadline, so the
    call is in flight for about a second. Meanwhile the main thread and a
    third thread pass an event back and forth. The handshake finishing
    before the computation stops proves the GIL was free."""
    ring = sylvester.PolynomialRing.prime_field(
        32003, [f"x{index}" for index in range(1, 9)]
    )
    ideal = ring.ideal(cyclic(ring, 8))

    running = threading.Event()
    stopped = threading.Event()
    ping = threading.Event()
    pong = threading.Event()
    outcome = []

    def compute():
        running.set()
        try:
            ideal.groebner_basis(timeout=1.0, threads=1)
            outcome.append("finished")
        except sylvester.BudgetExhausted:
            outcome.append("stopped")
        stopped.set()

    def answer():
        ping.wait()
        pong.set()

    responder = threading.Thread(target=answer)
    worker = threading.Thread(target=compute)
    responder.start()
    worker.start()

    running.wait()
    ping.set()
    assert pong.wait(timeout=60.0), "the handshake deadlocked behind the GIL"
    overlapped = not stopped.is_set()

    worker.join()
    responder.join()
    assert outcome == ["stopped"]
    assert overlapped


def test_two_threads_compute_at_once(cyclic):
    ring = sylvester.PolynomialRing.prime_field(32003, ["x", "y", "z"])
    ideal = ring.ideal(cyclic(ring, 3))
    expected = [str(f) for f in ideal.groebner_basis()]
    results = []
    lock = threading.Lock()

    def compute():
        basis = ideal.groebner_basis()
        with lock:
            results.append([str(f) for f in basis])

    threads = [threading.Thread(target=compute) for _ in range(4)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert results == [expected] * 4
