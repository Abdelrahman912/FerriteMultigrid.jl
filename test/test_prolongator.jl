function create_1d_1element_mass_matrix(p = 1, nqp = 2)
    grid = generate_grid(Line, (1,), Vec((0.0,)), Vec((1.0,)))

    ip = Lagrange{RefLine,p}()
    qr = QuadratureRule{RefLine}(nqp)
    cellvalues = CellValues(qr, ip)

    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    n_basefuncs = getnbasefunctions(cellvalues)
    Me = zeros(n_basefuncs, n_basefuncs)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        _element_mass_matrix!(Me, cellvalues)
    end
    return Me
end

@testset "Mass Matrix" begin
    # p =1 ,nqp = 2
    Me = create_1d_1element_mass_matrix(1, 2)
    Me_expected = 1/6 * [2 1; 1 2]
    @test Me ≈ Me_expected
    # p = 2, nqp = 3
    Me2 = create_1d_1element_mass_matrix(2, 3)
    Me2_expected = (1/30) * [
        4 -1 2;
        -1 4 2;
        2 2 16
    ]

    @test Me2 ≈ Me2_expected
end
