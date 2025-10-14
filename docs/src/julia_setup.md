# Julia Setup for PowerDynamics

In this document we'll provide a very brief overview on how to set up your development environment with Julia.
The goal of this document is to get you ready to follow our examples in this documentation.

We recommend you checking out further resources, like the excellent
- community provided [Modern Julia Workflows blog](https://modernjuliaworkflows.org/) and the
- [official Julia documentation](https://docs.julialang.org).

## Setup Julia
First you need to install Julia.
Check out the [install instructions](https://julialang.org/install/) on the official Julia homepage.
It is recommended installing Julia using **juliaup**.

!!! note
[juliaup](https://github.com/JuliaLang/juliaup) is a Julia version multiplexer, you can use it to "manage" different Julia version on the same system
using commands like `juliaup update` (update installed versions), `juliaup add 1.11` (add a new version), `juliaup default 1.11` (tell you computer which version to start when you invoke the `julia` command).

In general, julia can be run without any administrator privileges.

**Windows***
You can install Julia on Windows using the Microsoft Store. Instead of searching manually, you can invoke this installation in the terminal.
Search for `powershell` in the start menu and execute
```bash
winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
```

**Linux/Mac**
Execute
```bash
curl -fsSL https://install.julialang.org | sh
```
in your terminal.

The installer will ask you some questions, you won't need to change any defaults.

**Verify installation***
After the installation, run the `julia` command in your terminal. You should get an output like this:
![image](../assets/juliasetup/julia-prompt.png)

!!! tip "Stick to Julia 1.11 for now"
    Per default, juliaup will install the latest stable version of julia. As of 10/2025, that is 1.12. However there are some minor problems in 1.12 which might make it less stable. For now, I'd advise you to stay at Julia 1.11. To do so, run
    ```bash
    juliaup add 1.11
    juliaup default 1.11
    ```
    after which your `julia` command should point to Julia 1.11 per default.

## The REPL
Julias main interface is the REPL (read-eval-print-loop). It is similar to `ipython` or the Matlab command window.
You can execute code in the REPL by typing it. Try `julia> println("hello world")<RETURN>` to run your first Julia command.
Exit the REPL by pressing `<ctrl>+d` or typing `exit()<RETURN>`.

!!! tip
    The REPL has lots of great features. Just to name a few:
    - hit `]` to enter package manager mode (see [Environment Basics](@ref env-basics) below)
    - hit `?` to enter help mode: search for functions or concepts to get documentation, for example, look up the docstring of `println`
    - scroll through history using arrow up and arrow down
    The REPL is a great to get to know julia.

## Setup VSCode with the Julia extension
The REPL is great, but its not great for actual development. For that we need an editor.
The correct editor for most people will be [Visual Studio Code](https://code.visualstudio.com/).
Please [download and install](https://code.visualstudio.com/Download) it.

Within VSCode, you need to install the [julia-extension](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia).
You should be able to click the "install" link in any browser, it'll open VSCode and install the extension.

## [Environment and Package Manager Basics](@id env-basics)
Julia has a built in package manager. It lets you install, update and manage dependencies across your projects.
Often, different projects require different dependencies and different versions.

!!! tip "Example"
    Lets say you have a research project using `PowerDynamics` at version `v1.0.0` and its working great.
    However, you want to start a new project, which requires a new feature introduced in version `v2.2.0`!
    Environments help you to set up two different folders, one of which uses `v1.0.0` and one uses `v2.2.0`.
    This approach helps you to try out new stuff without the fear of breaking your already working projects
    due to incompatible updates.

To create your first environment, create a new folder somewhere (for example `/home/Documents/powerdynamics_playground`).
You can open this (empty) folder in VSCode and create a new file for your first script:
![image](../assets/juliasetup/vscode-script.png)

Open the newly created file, hit `<Ctrl> + <Shift> + P` to bring up the VSCode "command palette", search for "Start Repl"
and launch your Julia REPL.

In the REPL, you can execute `pwd()` (print working directory) to see the directory where your REPL was launched.
Hit `]` to launch the package manager. Your REPL changes:
```julia-repl
(@v1.11) pkg>
```
This output means, your "active" environment is the global environment for julia `v1.11`.
If you were to add packages here, you would add them globally. Normally, that's not what we want.
So instead we activate the current folder (`.`) as our working environment:
```julia-repl
(@v1.11) pkg> activate .
```
After activation, we can add `PowerDynamics` to our newly created environment using
```julia-repl
(powerdynamics_playground) pkg> add PowerDynamics
```
This will install PowerDynamics and all its dependencies. It will also precompile all of that. This may take a while...

![image](../assets/juliasetup/vscode-add-pd.png)

Once you've added a package, you'll see that two new files appeared: `Project.toml` and `Manifest.toml`
The `Project.toml` lists all of your dependencies. This is the file you change when "adding" new packages.
It is relatively short and human readable. It only lists the top level dependencies, in this case just `PowerDynamics`.
The `Manifest.toml` on the other hand lists *all* packages in the current environment. It is a full snapshot containing
all versions of every dependency and also the dependencies of your dependencies. Never edit it by hand!
```
ProjectRoot
├╴example_script.jl
├╴Project.toml
╰╴Manifest.toml
```
The existence of a `Project.toml` file marks a folder as a "project". You can "activate" such a folder.
When executing code in a file using the julia VSCode extension, it will automatically activate the environment for you.

## Executing Code in VSCode
Working with Julia is much like working in a notebook, due to the persistent REPL.
Because of how Julia works internally, everything you do will take much more time the first time you do it.
Therefore, it is always preferred to have a **persistent REPL** over relaunching julia. That is:
- opening a REPL,
- executing a script in REPL,
- changing script and
- run script again in same REPL
is much preferred to let's say
```bash
julia myscript.jl
```
where you would pay the startup costs on every iteration.

Besides executing an entire script (play button up top), you can execute **single lines** and **code blocks** in the VSCode:
- hover your cursor on a line and hit `<Shift>+<RETURN>` to "send" that line of code to the repl
- select multiple lines and hit `<SHIFT>+<RETURN>` to "send" all selected lines or
- hit `<ALT>+<RETURN>` (or `<CMD>+<RETURN>`) to execute entire *code cell*, where a code cell is everything between lines starting with `##`.

!!! note "Task"
    Install `ModelingToolkit` in addition to `PowerDynamics` in your environment, copy the following code to your script:
    ```julia
    using PowerDynamics
    using ModelingToolkit
    using PowerDynamics: Library

    @named swing = Library.Swing(V=1)
    busmodel = MTKBus(swing)
    swingbus = compile_bus(busmodel)
    ```
    execute it line by line and enjoy your first busmodel!
    ![image](../assets/juliasetup/vscode-firstmodel.png)
    
## Running a PowerDynamics Examples
Now you know everything you need to know to run our examples locally.

At the beginning of each example, there is a link:

> This tutorial can be downloaded as a normal Julia script here

Go to the [getting started](@ref getting-started) tutorial, download the script, put it in your
directory and go through it, executing it block by block.

!!! tip
    You can install multiple packages at once using commas:
    ```julia-repl
    julia> ] add OrdinaryDiffEqRosenbrock, CairoMakie
    ```

## Defining Functions and Revise.jl
Eventually, you'll "grow out" of putting everything into a script.
The simplemost next step is to put parts of you code in functions.
You could put those functions in the script, however if you change the code you need to evaluate them again.
Alternatively, you can create a new file, for example `myfunctions.jl`:
```julia
# contents of myfunctions.jl
function foo()
    println("Hello World")
end
```
Then, you can use [Revise.jl](https://github.com/timholy/Revise.jl) to track changes in that file.
In you main script write
```julia
Revise.includet("myfunctions.jl") # <- execute this line

foo() # <- prints "hello world"
```
then, you can update your `myfunctions.jl` file
```julia
# contents of myfunctions.jl
function foo()
    println("Hallo Welt")
end
```
and save it. If you evaluate `foo()` again (either in REPL directly or in your script), it'll print "Hallo Welt" -- Revise automatically updates the function definition.

At some point, you might want to put your project in a "Julia Package", but that's beyond the scope of this tutorial. Consult the documentation and your favorite Chatbot for help. Especially check out `] dev MyPackage` in contrast to `] add MyPackage`.

## Additional Remarks & Tips
1. Embrace the *interactivity of the REPL*: it really is (one of) Julias superpower. Play around with objects, make use of the help mode and use functions like `propertynames()`.
2. Get to know the package manager. It is really good! However, package management is not a trivial problem so package managers tend to be non-trivial. Read the docs and try to understand it. Just using it by trial and error will lead to major frustration down the road. Especially the option to `] dev` packages and define local [sources](https://pkgdocs.julialang.org/v1.12/toml-files/#The-%5Bsources%5D-section) can be very powerful tools.
3. Seek help. The [julia discourse](https://discourse.julialang.org/) and the [Julia Slack](https://julialang.org/slack/) and the [Julia Zulip](https://julialang.zulipchat.com/) are great platforms to post your problems and get help from the community.
