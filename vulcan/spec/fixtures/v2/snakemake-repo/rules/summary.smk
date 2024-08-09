configfile: "config.yaml"

rule summary:
    input:
        poem_count="output/count_poem.txt",
        poem_count_2="output/count_poem_2.txt",
        arithmetic="output/arithmetic.txt",
        checker="output/check.txt"
    output:
        "output/summary.txt"
    params:
        final_message=config["final_message"]
    run:
        with open(output[0], "w") as f:
            f.write(f"Count of poem.txt: {open(input.poem_count).read().strip()}\n")
            f.write(f"Count of poem_2.txt: {open(input.poem_count_2).read().strip()}\n")
            f.write(f"Result of arithmetic operation: {open(input.arithmetic).read().strip()}\n")
            f.write(f"Result of checker: {open(input.checker).read().strip()}\n")
            f.write(f"\n{params.final_message}\n")

