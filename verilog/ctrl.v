`timescale 1ns / 1ps

`define MEM_WIDTH 32
`define MEM_SIZE_WIDTH 11
`define MEM_SIZE (1<<`MEM_SIZE_WIDTH)

`define INS_WIDTH 16
`define PC_WIDTH (`MEM_SIZE_WIDTH/`INS_WIDTH)

`define GPR_SEL_WIDTH 3
`define WORD_WIDTH 256
`define GPR_SIZE (1<<`GPR_SEL_WIDTH)

//WRITEBACK info:
//write_double_word
//source[WB_SRC_WIDTH]
//dest[GPR_SEL_WIDTH]
`define WB_SRC_WIDTH 3
`define INS_WB_WIDTH (1+1+`WB_SRC_WIDTH+`GPR_SEL_WIDTH)
`define INS_EX_WIDTH (`INS_WB_WIDTH)

`define WB_SRC_NONE 0
`define WB_SRC_MEM 1
`define WB_SRC_ALU 2

`define ISA_OP_CLASS_WIDTH 2

`define ISA_OP_MEM_CLASS 0
`define ISA_OP_MEM_WIDTH 1
`define ISA_OP_MEM_R_WIDTH (`GPR_SEL_WIDTH)
`define ISA_OP_MEM_ADDR_WIDTH (`INS_WIDTH-`ISA_OP_CLASS_WIDTH-`ISA_OP_MEM_WIDTH-`ISA_OP_MEM_R_WIDTH-`GPR_SEL_WIDTH)

`define ISA_OP_ARITH_CLASS 1
`define ISA_OP_ARITH_WIDTH 3
`define ISA_OP_ARITH_RS0_WIDTH (`GPR_SEL_WIDTH)
`define ISA_OP_ARITH_RS1_WIDTH (`GPR_SEL_WIDTH)
`define ISA_OP_ARITH_RD0_WIDTH (`GPR_SEL_WIDTH)

`define ISA_OP_CTRL_CLASS 2
`define ISA_OP_CTRL_WIDTH 2
`define ISA_OP_CTRL_ADDR_WIDTH (`INS_WIDTH-`ISA_OP_CLASS_WIDTH-`ISA_OP_CTRL_WIDTH)

/*
module ctrl(
    input wire clk,
    input wire reset,
    input wire [`INS_WIDTH-1:0] ins_in,
    input wire [`MEM_WIDTH-1:0] dat_in,
    output reg [`MEM_WIDTH-1:0] dat_out,
    output reg [`ISA_OP_MEM_ADDR_WIDTH-1:0] dat_addr,
    output reg rd,
    output reg wr,
    output reg [`ISA_OP_CTRL_ADDR_WIDTH-1:0] ins_addr,
    output reg done
    );


reg [`PC_WIDTH-1:0] pc;
reg [`PC_WIDTH-1:0] next_pc;
reg fetch_stall;
reg ex_stall;
reg wb_stall;
reg [`INS_WIDTH-1:0] ins_fetch;
reg [`INS_EX_WIDTH-1:0] ins_ex;
reg [`INS_WB_WIDTH-1:0] ins_wb;
reg [`WORD_WIDTH-1:0] sbus;
wire [`WORD_WIDTH-1:0] dbus;

reg ins_ex_done;
reg do_wb;
reg ins_wb_done;

reg alu_add;
reg alu_shl;
reg alu_sto;
wire alu_eq;
wire alu_mz;
wire [`WORD_WIDTH-1:0] alu_out;
//reg mul_start;
//wire mul_done;
//wire [`WORD_WIDTH-1:0] mul_out;

always @* ins_fetch = ins_in;
always @* begin
    case(ins_fetch[0+:`ISA_OP_CLASS_WIDTH])
    `ISA_OP_MEM_CLASS: begin
    end
    `ISA_OP_ARITH_CLASS: begin
    end
    `ISA_OP_CTRL_CLASS: begin
    end
    endcase
end

always @* begin
    ins_ex_done = 1'b1;
    do_wb = 1'b0;
    next_pc = pc;
    case(ins_ex[0+:`WB_SRC_WIDTH])
    `WB_SRC_NONE: begin
    next_pc = pc + 1'b1;
    end
    `WB_SRC_MEM: begin
    end
    `WB_SRC_ALU: begin
    next_pc = pc + 1'b1;
    end
    endcase
end

always @* begin
    ins_wb_done = 1'b0;
    case(ins_wb[0+:`WB_SRC_WIDTH])
    `WB_SRC_NONE: begin
    ins_wb_done = 1'b1;
    end
    `WB_SRC_MEM: begin
    //ins_wb_done = mem_write_done;
    end
    `WB_SRC_ALU: begin
    ins_wb_done = 1;
    end
    endcase
end





always @(posedge clk, posedge reset) begin
    if(reset) begin
        pc <= 0;
    end else begin
        if(ins_ex_done) begin
            pc <= next_pc;
        end
    end
end



endmodule
*/

module data_path(
    input wire clk,
    input wire reset,
    input wire start,
    //operations using sbus (read from GPR)
    input wire alu_add,
    input wire alu_shl,
    input wire alu_sto,
    input wire mul_stos,
    input wire mul_stox,
    input wire div_stoy,
    input wire div_stox,
    input wire mem_wr,
    //operations using dbus (write to GPR)
    input wire alu_to_gpr,
    input wire mul_to_gpr,
    input wire div_to_gpr,
    input wire mem_to_gpr,
    //misc operations
    input wire mem_ld_poly,
    input wire mem_ld_poly_msb,
    input wire mul_mod,
    input wire div_mod,
    input wire mul_plain,
    input wire mul_clear,
    //operations parameters
    input wire mul_dbus_sel,
    input wire [`MEM_WIDTH-1:0] dat_in,
    input wire [`ISA_OP_MEM_ADDR_WIDTH-1:0] dat_addr_in,
    input wire [`GPR_SEL_WIDTH-1:0] gpr_rd_addr,
    input wire [`GPR_SEL_WIDTH-1:0] gpr_wr_addr,
    //outputs
    output wire alu_eq,
    output wire alu_early_mz,
    output wire alu_mz,
    output reg [`MEM_WIDTH-1:0] dat_out,
    output reg [`ISA_OP_MEM_ADDR_WIDTH-1:0] dat_addr_out,
    output reg wr,
    output reg run,
    output reg early_done,
    output reg done
    );

reg [`WORD_WIDTH-1:0] gpr[`GPR_SIZE:0];
wire [`WORD_WIDTH-1:0] sbus = gpr[gpr_rd_addr];
reg [`WORD_WIDTH-1:0] dbus;

reg [`WORD_WIDTH-1:0] irreducible_poly;
reg [`WORD_WIDTH-1:0] irreducible_poly_msb;

wire alu_done;
wire alu_run;
wire [`WORD_WIDTH-1:0] alu_out;

wire mul_run;
wire mul_done;
wire [`WORD_WIDTH-1:0] mul_out;

wire div_run;
wire div_done;
wire [`WORD_WIDTH-1:0] div_out;

reg gpr_wr_done;

reg mem_run;
reg mem_done;
//reg mem_run;
//wire mem_early_done = mem_run & (step == ((`WORD_WIDTH / `MEM_WIDTH)-1));
//reg mem_done;
//always @(posedge clk) begin
//    if(reset) begin
//        mem_run <= 0;
//        mem_done <= 0;
//    end else begin
//        if(start) mem_run <= mem_wr|mem_ld_poly|mem_ld_poly_msb|mem_to_gpr;
//        if(mem_early_done) mem_run <= 0;
//        mem_done <= mem_early_done;
//    end
//end
//always @* early_done = mem_early_done;

reg [7:0] step;

always @* begin
    mem_run = mem_wr|mem_ld_poly|mem_ld_poly_msb|mem_to_gpr;
    mem_done = mem_run & (step == (`WORD_WIDTH / `MEM_WIDTH));
    run = div_run | mul_run | alu_run | mem_run;
    done = div_done | mul_done | alu_done | gpr_wr_done | mem_done;
end

always @* dat_addr_out = dat_addr_in + step;
always @* dat_out = sbus[step*`MEM_WIDTH+:`MEM_WIDTH];
always @* wr = mem_wr;

always @* begin
    gpr_wr_done = 0;
    case(1'b1)
    alu_to_gpr: gpr_wr_done = 1'b1;
    mul_to_gpr: gpr_wr_done = 1'b1;
    div_to_gpr: gpr_wr_done = 1'b1;
    mem_to_gpr: gpr_wr_done = step == (`WORD_WIDTH / `MEM_WIDTH);
    endcase
end

always @(posedge clk) begin: STEP
    if(reset) step <= 0;
    else begin
        if(run) step <= step + 1;
        if(done) step <= 0;
    end
end

always @(posedge clk) begin: GPR
    case(1'b1)
    alu_to_gpr: gpr[gpr_wr_addr] <= alu_out;
    mul_to_gpr: gpr[gpr_wr_addr] <= mul_out;
    div_to_gpr: gpr[gpr_wr_addr] <= div_out;
    mem_to_gpr: gpr[gpr_wr_addr][step*`MEM_WIDTH+:`MEM_WIDTH] <= dat_in;
    endcase
end

always @(posedge clk) begin: POLY
    if(reset)
        irreducible_poly <= 0;
    else begin
        if(mem_ld_poly) begin
            irreducible_poly[step*`MEM_WIDTH+:`MEM_WIDTH] <= dat_in;
        end
    end
end

always @(posedge clk) begin: POLY_MSB
    if(reset)
        irreducible_poly_msb <= 0;
    else begin
        if(mem_ld_poly_msb) begin
            irreducible_poly_msb[step*`MEM_WIDTH+:`MEM_WIDTH] <= dat_in;
        end
    end
end

gf2m_alu u_g2m_alu(
    .clk(clk),.reset(reset),
    .add(alu_add),.shl(alu_shl),.sto(alu_sto),
    .sbus(sbus),.dbus(alu_out),
    .eq(alu_eq),.early_mz(alu_early_mz),.mz(alu_mz),
    .run(alu_run),.done(alu_done));
gf2m_mul u_g2m_mul(
    .clk(clk),.reset(reset),
    .mod_mul(mul_mod),
    .plain_mul(mul_plain),
    .stos(mul_stos),
    .stox(mul_stox),
    .dbus_sel(mul_dbus_sel),
    .clear(mul_clear),
    .irreducible_poly_msb(irreducible_poly_msb),
    .irreducible_poly(irreducible_poly),
    .sbus(sbus),.dbus(mul_out),
    .run(mul_run),.done(mul_done)
);
gf2m_div u_gf2m_div(
    .clk(clk),.reset(reset),
    .mod_div(div_mod),//could also support inversion almost for free
    .stoy(div_stoy),//dividend
    .stox(div_stox),//divider
    .irreducible_poly(irreducible_poly),
    .sbus(sbus),.dbus(div_out),
    .run(div_run),.done(div_done)
    );
endmodule
