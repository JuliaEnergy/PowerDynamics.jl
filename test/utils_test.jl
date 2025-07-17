using PowerDynamics
using PowerDynamics.Library
using PowerDynamics: pin_parameters
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using ModelingToolkit: get_name

@info "Start Utils tests"

@testset "test @attach_metadata!" begin
    @mtkmodel MyModel begin
        @components begin
            busbar = BusBar()
        end
    end
    @test_throws UndefKeywordError MyModel()
    @test get_name(MyModel(name=:bar)) == :bar
    @test ModelingToolkit.get_metadata(MyModel(name=:bar)) == nothing
    # @b MyModel(name=:bar)

    @attach_metadata! MyModel (;origin="file", name=:defname)
    @test ModelingToolkit.get_metadata(MyModel()) == (;origin="file")
    @test get_name(MyModel()) == :defname
    @test get_name(MyModel(name=:bar)) == :bar
    # @b MyModel()

    @attach_metadata! MyModel (;another=:foo)
    @test ModelingToolkit.get_metadata(MyModel()) == (;origin="file", another=:foo)
    @test get_name(MyModel()) == :defname

    @attach_metadata! MyModel (;name=:newname)
    @test ModelingToolkit.get_metadata(MyModel()) == (;origin="file", another=:foo)
    @test get_name(MyModel()) == :newname
end

@testset "pin parameters" begin
    @mtkmodel Inner begin
        @parameters begin
            a = 0
            b = 0
        end
        @variables begin
            i(t)
        end
        @equations begin
            i ~ a + b
        end
    end
    @mtkmodel Outer begin
        @components begin
            inner = Inner()
        end
        @parameters begin
            a = 0
            b = 0
        end
        @variables begin
            o(t)
        end
        @equations begin
            o ~ inner.i * (a + b)
        end
    end
    @named outer = Outer()
    full_equations(outer)
    sys1 = pin_parameters(outer, :a => 1)
    sys2 = pin_parameters(outer, :outer₊b => 1)
    full_equations(sys1)
    full_equations(sys2)
    sys1 = pin_parameters(outer, outer.inner.a => 1)
    sys2 = pin_parameters(outer, :inner₊a => 1)
    @test full_equations(sys1) == full_equations(sys2)
end
