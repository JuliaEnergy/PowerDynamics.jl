using PowerDynamics
using PowerDynamics.Library
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using ModelingToolkit: get_name
using PowerDynamics: get_custom_metadata

@info "Start Utils tests"

@testset "test @attach_metadata!" begin
    @mtkmodel MyModel begin
        @components begin
            busbar = BusBar()
        end
    end
    @test_throws UndefKeywordError MyModel()
    @test get_name(MyModel(name=:bar)) == :bar
    default_metadata = ModelingToolkit.get_metadata(MyModel(name=:bar))
    # @b MyModel(name=:bar)

    @attach_metadata! MyModel (;origin="file", name=:defname)
    @test get_custom_metadata(MyModel(), :origin) == "file"
    @test get_name(MyModel()) == :defname
    @test get_name(MyModel(name=:bar)) == :bar
    # @b MyModel()

    @attach_metadata! MyModel (;another=:foo)
    @test get_custom_metadata(MyModel()) == (; origin = "file", another = :foo)
    @test get_name(MyModel()) == :defname

    @attach_metadata! MyModel (;name=:newname)
    @test get_custom_metadata(MyModel()) == (; origin = "file", another = :foo)
    @test get_name(MyModel()) == :newname
end
