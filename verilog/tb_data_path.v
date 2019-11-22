`timescale 1ns / 1ps


`define MEM_WIDTH 32
//`define MEM_SIZE_WIDTH 11
//`define MEM_SIZE (1<<`MEM_SIZE_WIDTH)

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

`define IMEM_WIDTH `INS_WIDTH
`define IMEM_ADDR_WIDTH `ISA_OP_CTRL_ADDR_WIDTH

`define DMEM_WIDTH `MEM_WIDTH
`define DMEM_ADDR_WIDTH `ISA_OP_MEM_ADDR_WIDTH

`define assert(signal, value) \
        if (signal !== value) begin \
            $display("%t ASSERTION FAILED in %m: signal != value",$time); \
            wait_clocks(10);\
            $finish; \
        end

//synchronous single port RAM
module ram#(
    parameter WIDTH = 32,
    parameter ADDR_WIDTH = 9
)(
    input wire clk,
    input wire wr,
    input wire [ADDR_WIDTH-1:0] addr,
    input wire [WIDTH-1:0] din,
    output reg [WIDTH-1:0] dout
);
localparam DEPTH = (1<<ADDR_WIDTH);
reg [WIDTH-1:0] storage[DEPTH-1:0];
always @(posedge clk) begin
    if(wr) storage[addr] <= din;
end
always @* dout = storage[addr];
endmodule



module tb_data_path();
wire dbg_print_gpr = 0;
localparam IMEM_WIDTH = `IMEM_WIDTH;
localparam IMEM_ADDR_WIDTH = `IMEM_ADDR_WIDTH;

localparam DMEM_WIDTH = `DMEM_WIDTH;
localparam DMEM_ADDR_WIDTH = `DMEM_ADDR_WIDTH;

localparam WORD_WIDTH = `WORD_WIDTH;
localparam DMEM_PER_WORD = WORD_WIDTH / DMEM_WIDTH;

reg clk;
reg reset;
reg dmem_wr;
reg [DMEM_ADDR_WIDTH-1:0]dmem_addr;
reg [DMEM_WIDTH-1:0]dmem_din;
wire [DMEM_WIDTH-1:0]dmem_dout;

wire imem_wr;
wire [IMEM_ADDR_WIDTH-1:0]imem_addr;
wire [IMEM_WIDTH-1:0]imem_din;
wire [IMEM_WIDTH-1:0]imem_dout;

ram #(.WIDTH(DMEM_WIDTH),.ADDR_WIDTH(DMEM_ADDR_WIDTH)) dmem(.clk(clk), .wr(dmem_wr), .addr(dmem_addr), .din(dmem_din), .dout(dmem_dout));
ram #(.WIDTH(IMEM_WIDTH),.ADDR_WIDTH(IMEM_ADDR_WIDTH)) imem(.clk(clk), .wr(imem_wr), .addr(imem_addr), .din(imem_din), .dout(imem_dout));

reg start;
reg alu_add;
reg alu_shl;
reg alu_sto;
reg mul_stos;
reg mul_stox;
reg div_stoy;
reg div_stox;
reg mem_wr;
//operations using dbus (write to GPR)
reg alu_to_gpr;
reg mul_to_gpr;
reg div_to_gpr;
reg mem_to_gpr;
//misc operations
reg mem_ld_poly;
reg mem_ld_poly_msb;
reg mul_mod;
reg mul_plain;
reg mul_clear;
reg div_mod;
//operations parameters
reg mul_dbus_sel;
reg [`MEM_WIDTH-1:0] dat_in;
reg [`ISA_OP_MEM_ADDR_WIDTH-1:0] dat_addr_in;
reg [`GPR_SEL_WIDTH-1:0] gpr_rd_addr;
reg [`GPR_SEL_WIDTH-1:0] gpr_wr_addr;
//outputs
wire alu_eq;
wire alu_early_mz;
wire alu_mz;
wire [`MEM_WIDTH-1:0] dat_out;
wire [`ISA_OP_MEM_ADDR_WIDTH-1:0] dat_addr_out;
wire wr;
wire run;
wire early_done;
wire done;

data_path u_data_path(
    .clk(clk),.reset(reset),.start(start),
    .alu_add(alu_add),
    .alu_shl(alu_shl),
    .alu_sto(alu_sto),
    .mul_stos(mul_stos),
    .mul_stox(mul_stox),
    .div_stoy(div_stoy),
    .div_stox(div_stox),
    .mem_wr(mem_wr),
    .alu_to_gpr(alu_to_gpr),
    .mul_to_gpr(mul_to_gpr),
    .div_to_gpr(div_to_gpr),
    .mem_to_gpr(mem_to_gpr),
    .mem_ld_poly(mem_ld_poly),
    .mem_ld_poly_msb(mem_ld_poly_msb),
    .mul_mod(mul_mod),
    .mul_plain(mul_plain),
    .mul_clear(mul_clear),
    .mul_dbus_sel(mul_dbus_sel),
    .div_mod(div_mod),
    .dat_in(dat_in),
    .dat_addr_in(dat_addr_in),
    .gpr_rd_addr(gpr_rd_addr),
    .gpr_wr_addr(gpr_wr_addr),
    .alu_eq(alu_eq),
    .alu_mz(alu_mz),
    .alu_early_mz(alu_early_mz),
    .dat_out(dat_out),
    .dat_addr_out(dat_addr_out),
    .wr(wr),
    .run(run),
    .early_done(early_done),
    .done(done)
);

always @* begin
    dmem_addr = dat_addr_out;
    dmem_din = dat_out;
    dat_in = dmem_dout;
    dmem_wr = wr;
end


task wait_clock;
begin
    @ (negedge clk);
end
endtask

task reset_dut;
begin
    wait_clock();
    reset = 1;
    wait_clock();
    reset = 0;
end
endtask

task wait_clocks;
    input integer nclocks;
    integer i;
begin
    for(i=0;i < nclocks; i= i+1) begin
        wait_clock();
    end
end
endtask

task start_op;
begin
    start = 1;
    @(negedge clk);
    start = 0;
end
endtask

reg done_reg;
always @(posedge clk) done_reg <= early_done;

task wait_done;
begin
    wait_clock();
    while(!done) begin
        wait_clock();
    end
    //while(!done_reg) begin
    //    @ (posedge clk);
    //end
end
endtask

task write_dmem;
    input [DMEM_ADDR_WIDTH-1:0] addr;
    input [DMEM_WIDTH-1:0] dat;
begin
    dmem.storage[addr] = dat;
end
endtask

task read_dmem;
    input [DMEM_ADDR_WIDTH-1:0] addr;
    output [DMEM_WIDTH-1:0] dat;
begin
    dat = dmem.storage[addr];
end
endtask

`define check_dmem(addr,dat) `assert(dmem.storage[addr],dat);
//task check_dmem;
//    input [DMEM_ADDR_WIDTH-1:0] addr;
//    input [DMEM_WIDTH-1:0] dat;
//begin
//    `assert(dmem.storage[addr],dat);
//end
//endtask

task check_alu_eq;
    input dat;
begin
    `assert(alu_eq,dat);
end
endtask

task check_alu_mz;
    input dat;
begin
    `assert(alu_mz,dat);
end
endtask

task check_gpr;
    input [`GPR_SEL_WIDTH-1:0] gpr;
    input [WORD_WIDTH-1:0] dat;
begin
    //$display(u_data_path.gpr[gpr]);
    //$display(dat);
    `assert(u_data_path.gpr[gpr],dat);
end
endtask

task disp_gpr;
    input [`GPR_SEL_WIDTH-1:0] gpr;
begin
    if(dbg_print_gpr) $display("gpr[%d]=0x%x",gpr,u_data_path.gpr[gpr]);
end
endtask

task check_poly;
    input [WORD_WIDTH-1:0] dat;
begin
    `assert(dat,u_data_path.irreducible_poly);
end
endtask

task check_poly_msb;
    input [WORD_WIDTH-1:0] dat;
begin
    `assert(dat,u_data_path.irreducible_poly_msb);
end
endtask

task write_dmem_word;
    input [DMEM_ADDR_WIDTH-1:0] addr;
    input [WORD_WIDTH-1:0] dat;
    reg [DMEM_ADDR_WIDTH-1:0] mem_addr;
    integer i;
begin
    mem_addr = addr * DMEM_PER_WORD;
    for(i=0;i<DMEM_PER_WORD;i=i+1) begin
        write_dmem(mem_addr+i,dat[i*DMEM_WIDTH+:DMEM_WIDTH]);
    end
end
endtask

task read_dmem_word;
    input [DMEM_ADDR_WIDTH-1:0] addr;
    output [WORD_WIDTH-1:0] dat;
    reg [DMEM_ADDR_WIDTH-1:0] mem_addr;
    integer i;
begin
    mem_addr = addr * DMEM_PER_WORD;
    for(i=0;i<DMEM_PER_WORD;i=i+1) begin
        read_dmem(mem_addr+i,dat[i*DMEM_WIDTH+:DMEM_WIDTH]);
    end
end
endtask

task check_dmem_word;
    input [DMEM_ADDR_WIDTH-1:0] addr;
    input [WORD_WIDTH-1:0] dat;
    reg [DMEM_ADDR_WIDTH-1:0] mem_addr;
    integer i;
begin
    mem_addr = addr * DMEM_PER_WORD;
    for(i=0;i<DMEM_PER_WORD;i=i+1) begin
        `check_dmem(mem_addr+i,dat[i*DMEM_WIDTH+:DMEM_WIDTH]);
    end
end
endtask

task ins_mem_ld_poly;
    input [DMEM_ADDR_WIDTH-1:0] addr;
begin
    dat_addr_in = addr * DMEM_PER_WORD;
    mem_ld_poly = 1;
    wait_done();
    mem_ld_poly = 0;
end
endtask

task ins_mem_ld_poly_msb;
    input [DMEM_ADDR_WIDTH-1:0] addr;
begin
    dat_addr_in = addr * DMEM_PER_WORD;
    mem_ld_poly_msb = 1;
    wait_done();
    mem_ld_poly_msb = 0;
end
endtask

task ins_mem_ld;
    input [DMEM_ADDR_WIDTH-1:0] addr;
    input [`GPR_SEL_WIDTH-1:0] dst;
begin
    @(negedge clk);
    dat_addr_in = addr * DMEM_PER_WORD;
    gpr_wr_addr = dst;
    mem_to_gpr = 1;
    //start_op();
    wait_done();
    mem_to_gpr = 0;
end
endtask

task ins_mem_st;
    input [DMEM_ADDR_WIDTH-1:0] addr;
    input [`GPR_SEL_WIDTH-1:0] src;
begin
    @(negedge clk);
    dat_addr_in = addr * DMEM_PER_WORD;
    gpr_rd_addr = src;
    mem_wr = 1;
    //start_op();
    wait_done();
    mem_wr = 0;
end
endtask

task ins_add;
    input [`GPR_SEL_WIDTH-1:0] dst;
    input [`GPR_SEL_WIDTH-1:0] a;
    input [`GPR_SEL_WIDTH-1:0] b;
begin
    gpr_rd_addr = a;
    alu_sto = 1;
    wait_done();
    alu_sto = 0;
    gpr_rd_addr = b;
    alu_add = 1;
    wait_done();
    alu_add = 0;
    gpr_wr_addr = dst;
    alu_to_gpr = 1;
    wait_done();
    alu_to_gpr = 0;
    $display("add %d, %d, %d",dst,a,b);
    disp_gpr(dst);
end
endtask

task ins_shl;
    input [`GPR_SEL_WIDTH-1:0] dst;
    input [`GPR_SEL_WIDTH-1:0] a;
begin
    gpr_rd_addr = a;
    alu_shl = 1;
    wait_done();
    alu_shl = 0;
    gpr_wr_addr = dst;
    alu_to_gpr = 1;
    wait_done();
    alu_to_gpr = 0;
    $display("shl %d, %d",dst,a);
    disp_gpr(dst);
end
endtask

task ins_test;//set alu_mz and alu_eq according to the source reg
    input [`GPR_SEL_WIDTH-1:0] a;
begin
    gpr_rd_addr = a;
    alu_sto = 1;
    wait_done();
    alu_sto = 0;
    $display("test %d",a);
    disp_gpr(a);
end
endtask

task ins_mov;
    input [`GPR_SEL_WIDTH-1:0] dst;
    input [`GPR_SEL_WIDTH-1:0] src;
begin
    gpr_rd_addr = src;
    alu_sto = 1;
    wait_done();
    alu_sto = 0;
    gpr_wr_addr = dst;
    alu_to_gpr = 1;
    wait_done();
    alu_to_gpr = 0;
end
endtask

task ins_mul_plain;
    input [`GPR_SEL_WIDTH-1:0] dst;
    input [`GPR_SEL_WIDTH-1:0] s;
    input [`GPR_SEL_WIDTH-1:0] x;
begin
    gpr_rd_addr = s;
    mul_stos = 1;
    wait_done();
    mul_stos = 0;
    gpr_rd_addr = x;
    mul_stox = 1;
    mul_clear = 1;
    wait_done();
    mul_stox = 0;
    mul_clear = 0;
    mul_plain = 1;
    wait_done();
    mul_plain = 0;
    gpr_wr_addr = dst;
    mul_dbus_sel=0;
    mul_to_gpr = 1;
    wait_done();
    mul_dbus_sel=1;
    gpr_wr_addr = dst+1;
    wait_done();
    mul_to_gpr = 0;
end
endtask

task ins_mul_mod;
    input [`GPR_SEL_WIDTH-1:0] dst;
    input [`GPR_SEL_WIDTH-1:0] s;
    input [`GPR_SEL_WIDTH-1:0] x;
begin
    gpr_rd_addr = s;
    mul_stos = 1;
    wait_done();
    mul_stos = 0;
    gpr_rd_addr = x;
    mul_stox = 1;
    mul_clear = 1;
    wait_done();
    mul_stox = 0;
    mul_clear = 0;
    mul_mod = 1;
    wait_done();
    mul_mod = 0;
    gpr_wr_addr = dst;
    mul_dbus_sel=0;
    mul_to_gpr = 1;
    wait_done();
    mul_to_gpr = 0;
    $display("mul_mod %d, %d, %d",dst,s,x);
    disp_gpr(dst);
end
endtask

task ins_div_mod;
    input [`GPR_SEL_WIDTH-1:0] dst;
    input [`GPR_SEL_WIDTH-1:0] y;
    input [`GPR_SEL_WIDTH-1:0] x;
begin
    gpr_rd_addr = y;
    div_stoy = 1;
    wait_done();
    div_stoy = 0;
    gpr_rd_addr = x;
    div_stox = 1;
    wait_done();
    div_stox = 0;
    div_mod = 1;
    wait_done();
    div_mod = 0;
    gpr_wr_addr = dst;
    div_to_gpr = 1;
    wait_done();
    div_to_gpr = 0;
    $display("div_mod %d, %d, %d",dst,y,x);
    disp_gpr(dst);
end
endtask


wire [234-1:0] modulus =         234'h20000000000000000000000000000000000000004000000000000000001;
wire [233-1:0] base_point_x    = 233'h17232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126;
wire [233-1:0] base_point_y    = 233'h1db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3;
wire [232-1:0] base_point_order = 232'h8000000000000000000000000000069d5bb915bcd46efb1ad5f173abdf;
wire [229-1:0] private_key = 229'h103b2142bdc2a3c3b55080d09df1808f79336da2399f5ca7171d1be9b0;
wire [231-1:0] public_key_x  = 231'h682886f36c68473c1a221720c2b12b9be13458ba907e1c4736595779f2;
wire [233-1:0] public_key_y = 233'h1b20639b41be0927090999b7817a3b3928d20503a39546044ec13a10309;
wire [231-1:0] k = 231'h73552f9cac5774f74f485fa253871f2109a0c86040552eaa67dba92dc9;
wire [232-1:0] kG_x = 232'hb8ad9c1d2cb29906e7d63c24601acbf4926d4eb7e9630458422469e5ef;
wire [230-1:0] r = 230'h38ad9c1d2cb29906e7d63c24601ac55736b438fb14f4093d6c32f63a10;
wire [231-1:0] s1 = 231'h7ca896d9c942e3246703e49de6e9e648b86493979e98f67bd1c1492705;
wire [230-1:0] s2 = 230'h2bd472bb73de51e649b1c6747bddf8d51b7585205ea4895b87ef9f59a9;
wire [231-1:0] s = 231'h647aad2599c21b6ee89be7ff957d98f684b7921de1fd3cc82c079624f4;

task load_modulus;
    input [DMEM_ADDR_WIDTH-1:0] dmem_addr;
    input [WORD_WIDTH-1:0] modulus;
    reg [WORD_WIDTH-1:0] tmp;
    integer i;
begin
    write_dmem_word(dmem_addr,modulus);
    ins_mem_ld_poly(dmem_addr);
    tmp = 1<<(WORD_WIDTH-1);
    while(0==(tmp & modulus)) tmp = tmp>>1;
    write_dmem_word(dmem_addr,tmp);
    ins_mem_ld_poly_msb(dmem_addr);
end
endtask

task init_curve_k233;
begin
    {start,alu_add,alu_shl,alu_sto,div_stoy,div_stox,div_mod,mul_stos,mul_stox,mul_mod,mul_plain,mul_clear,mul_dbus_sel,alu_to_gpr,div_to_gpr,mul_to_gpr,mem_to_gpr,mem_wr,mem_ld_poly,mem_ld_poly_msb} = 0;
    wait_clocks(10);
    curve_a2 = 0;
    write_dmem_word(0,modulus);
    write_dmem_word(1,base_point_x);
    write_dmem_word(2,base_point_y);
    write_dmem_word(3,base_point_order);
    write_dmem_word(4,private_key);
    write_dmem_word(5,k);
    //write_dmem_word(6,modulus & 256'hFFFFFFFF00000000000000000000000000000000000000000000000000000000);
    reset_dut();
    load_modulus(6,modulus);
    wait_clocks(10);
end
endtask

task test_mem_ops;
begin
    init_curve_k233();

    ins_mem_ld(1,0);
    ins_mem_ld(2,1);
    ins_mem_ld(4,2);

    ins_mem_ld_poly(0);
    ins_mem_ld_poly_msb(6);
    //`assert(modulus,u_data_path.irreducible_poly);
    //`assert(u_data_path.irreducible_poly_msb,256'h20000000000000000000000000000000000000000000000000000000000);
    check_poly(modulus);
    check_poly_msb(256'h20000000000000000000000000000000000000000000000000000000000);



    ins_mem_st(6,0);
    ins_mem_st(7,1);
    ins_mem_st(8,2);

    check_dmem_word(0,modulus);
    check_dmem_word(1,base_point_x);
    check_dmem_word(2,base_point_y);
    check_dmem_word(3,base_point_order);
    check_dmem_word(4,private_key);
    check_dmem_word(5,k);

    check_dmem_word(6,base_point_x);
    check_dmem_word(7,base_point_y);
    check_dmem_word(8,private_key);
    $display("test_mem_ops pass");
    wait_clocks(10);
end
endtask

task test_alu_add;
begin
    init_curve_k233();
    ins_mem_ld(0,0);
    check_gpr(0,modulus);
    ins_add(0,0,0);
    check_gpr(0,0);
    check_alu_eq(1);
    ins_mem_ld(0,1);
    ins_add(0,0,1);
    check_gpr(1,modulus);
    check_alu_eq(0);
    ins_mem_ld(1,1);
    ins_mem_ld(2,2);
    ins_add(7,2,1);
    check_gpr(7,base_point_x^base_point_y);
    $display("test_alu_add pass");
end
endtask

task test_alu_shl;
    input [`GPR_SEL_WIDTH-1:0] test_gpr;
begin
    init_curve_k233();
    ins_mem_ld(2,test_gpr);
    check_gpr(test_gpr,base_point_y);
    ins_shl(test_gpr,test_gpr);
    check_gpr(test_gpr,base_point_y<<1);
    check_alu_eq(0);
    check_alu_mz(0);
    ins_shl(test_gpr^1,test_gpr);
    check_gpr(test_gpr^1,base_point_y<<2);
    check_alu_eq(0);
    check_alu_mz(0);
    write_dmem_word(0,3<<(WORD_WIDTH-3));
    ins_mem_ld(0,test_gpr);
    ins_shl(test_gpr,test_gpr);
    check_alu_eq(0);
    check_alu_mz(0);
    ins_shl(test_gpr,test_gpr);
    check_alu_eq(0);
    check_alu_mz(1);
    ins_shl(test_gpr,test_gpr);
    check_alu_eq(1);
    check_alu_mz(1);
    ins_shl(test_gpr,test_gpr);
    check_alu_eq(1);
    check_alu_mz(0);
    $display("test_alu_shl pass");
end
endtask

task test_mul_plain_core;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_o;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_s;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_x;
    input [WORD_WIDTH-1:0] s;
    input [WORD_WIDTH-1:0] x;
    input [2*WORD_WIDTH-1:0] expected;
begin
    write_dmem_word(0,s);
    write_dmem_word(1,x);
    ins_mem_ld(0,test_gpr_s);
    ins_mem_ld(1,test_gpr_x);
    ins_mul_plain(test_gpr_o,test_gpr_s,test_gpr_x);
    check_gpr(test_gpr_o,expected[0+:WORD_WIDTH]);
    check_gpr(test_gpr_o+1,expected[WORD_WIDTH+:WORD_WIDTH]);
    ins_mem_ld(1,test_gpr_s);
    ins_mem_ld(0,test_gpr_x);
    ins_mul_plain(test_gpr_o,test_gpr_s,test_gpr_x);
    check_gpr(test_gpr_o,expected[0+:WORD_WIDTH]);
    check_gpr(test_gpr_o+1,expected[WORD_WIDTH+:WORD_WIDTH]);
end
endtask

task test_mul_mod_core;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_o;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_s;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_x;
    input [WORD_WIDTH-1:0] s;
    input [WORD_WIDTH-1:0] x;
    input [WORD_WIDTH-1:0] expected;
begin
    write_dmem_word(0,s);
    write_dmem_word(1,x);
    ins_mem_ld(0,test_gpr_s);
    ins_mem_ld(1,test_gpr_x);
    ins_mul_mod(test_gpr_o,test_gpr_s,test_gpr_x);
    check_gpr(test_gpr_o,expected);
    ins_mem_ld(1,test_gpr_s);
    ins_mem_ld(0,test_gpr_x);
    ins_mul_mod(test_gpr_o,test_gpr_s,test_gpr_x);
    check_gpr(test_gpr_o,expected);
end
endtask

task test_div_mod_core;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_o;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_y;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_x;
    input [WORD_WIDTH-1:0] y;
    input [WORD_WIDTH-1:0] x;
    input [WORD_WIDTH-1:0] expected;
begin
    write_dmem_word(0,y);
    write_dmem_word(1,x);
    ins_mem_ld(0,test_gpr_y);
    ins_mem_ld(1,test_gpr_x);
    ins_div_mod(test_gpr_o,test_gpr_y,test_gpr_x);
    check_gpr(test_gpr_o,expected);
end
endtask

task test_mul_plain;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_o;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_s;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_x;
begin
    test_mul_plain_core(test_gpr_o,test_gpr_s,test_gpr_x,4'h8,8'he7,11'h738);
    test_mul_plain_core(test_gpr_o,test_gpr_s,test_gpr_x,2'h3,14'h3b25,15'h4d6f);
    test_mul_plain_core(test_gpr_o,test_gpr_s,test_gpr_x,3'h6,230'h3f71bc32965b4dabba50d5a71d530eaed98127689d21663ab599cbcce9,232'h832588af75dbadf99de2fdd24fea27e6d506d3734ec7549fbd54b8aa76);
    test_mul_plain_core(test_gpr_o,test_gpr_s,test_gpr_x,256'h17232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126,256'h1db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3,512'h188cbf3002f0fceb9905f445cefa0c661ddaf85088e234de26a723857e03aa103f0259eaf60e30a4cf2404f9d08ceea4ad7a9200e8f6eab0f80aa);
    test_mul_plain_core(test_gpr_o,test_gpr_s,test_gpr_x,256'hF2C295F27A96B9435935807A7359F67FD014F9A8C9EE2589E13F0CC8B6630CA6,256'h876E46A6F24CE78C4D904AD897ECC395A14F3DFE78E803FC10D5A8DF4C632923,512'h7b9a69879baeac0dc1c9bb962d704bfaa0b5919ab00d1a2227b397211e40747634978dc6bc0cbde09bc5b93a7a8594f15272fcaf1aad7238b8fbcc5d8c93d72a);
    test_mul_plain_core(test_gpr_o,test_gpr_s,test_gpr_x,256'hFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,256'hFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,512'h55555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555);

    $display("test_mul_plain pass");
end
endtask

task test_mul_mod;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_o;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_s;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_x;
begin
    load_modulus(6,256'h11B);
    test_mul_mod_core(test_gpr_o,test_gpr_s,test_gpr_x,4'h2,8'h87,8'h15);
    load_modulus(6,modulus);
    test_mul_mod_core(test_gpr_o,test_gpr_s,test_gpr_x,256'h17232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126,256'h1db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3,256'h404c43af73958b87742ff9e35ec83a50fb77c1d266fa5b7e749ddd12ca);
    $display("test_mul_mod pass");
end
endtask

task test_div_mod;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_o;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_y;
    input [`GPR_SEL_WIDTH-1:0] test_gpr_x;
begin
    load_modulus(6,256'h13);
    test_div_mod_core(test_gpr_o,test_gpr_y,test_gpr_x,8'h8,8'h6,8'hd);
    test_div_mod_core(test_gpr_o,test_gpr_y,test_gpr_x,8'h5,8'hc,8'h4);
    load_modulus(6,modulus);
    test_div_mod_core(test_gpr_o,test_gpr_y,test_gpr_x,256'h183a4b8cc1296a55991b576cec7044c7dcf06dcfdef31bc10a569a1b53c,256'h157441042a78f77e1b884a309da0db7f3c42524d235e5431ab4b0acc123,256'h1b3f072ea66d0c24f08ae1f7b862af156be7542016a669c68c370e2fc73);
    test_div_mod_core(test_gpr_o,test_gpr_y,test_gpr_x,256'h1d5dcd80763c2fbc63292447884f7cdd7070395c41fb4202961599bf9a7,256'h5cecc597155e72ff4c971c270bffe9553acea85dbf515d1d01ce8c38a1,256'h1ebdb5bd60ceffecee2dff496b4dd7340e717796da8fb2ccb20229c30d4);
    $display("test_div_mod pass");
end
endtask

wire [WORD_WIDTH-1:0] origin = 1<<(WORD_WIDTH-1);

reg curve_a2;

task is_origin;
    input [`GPR_SEL_WIDTH-1:0] a;
    output o;
begin
    ins_test(a);
    o = alu_mz;
end
endtask

//point in r0,r1
//result in r0,r1
task point_double;
    reg flag;
    reg [`GPR_SEL_WIDTH-1:0] px;
    reg [`GPR_SEL_WIDTH-1:0] py;
    reg [`GPR_SEL_WIDTH-1:0] l;
begin
    px=0;//point.x in r0
    py=1;//point.y in r1
    l=2;
    is_origin(px,flag);
    if(!flag) begin
        ins_div_mod(l,py,px);
        ins_mul_mod(py,px,px);
        ins_add(l,l,px);
        ins_mul_mod(px,l,l);
        ins_add(px,px,l);
        if(curve_a2) begin
            //ins_add(px,px)//TODO: ins_inc
            `assert(1,0);
        end
        ins_mul_mod(l,l,px);
        ins_add(py,py,l);
        ins_add(py,py,px);
    end
end
endtask

//point a in r0,r1
//point b in r2,r3
//result in r0,r1
task point_add;
    reg flag;
    reg [`GPR_SEL_WIDTH-1:0] ax;
    reg [`GPR_SEL_WIDTH-1:0] ay;
    reg [`GPR_SEL_WIDTH-1:0] bx;
    reg [`GPR_SEL_WIDTH-1:0] by;
    reg [`GPR_SEL_WIDTH-1:0] l;
begin
    ax=0;//a.x in r0
    ay=1;//a.y in r1
    bx=2;//b.x in r2
    by=3;//b.y in r3
    l=4;
    /*is_origin(px,flag);
    if(!flag) begin
        ins_div_mod(l,py,px);
        ins_mul_mod(py,px,px);
        ins_add(l,l,px);
        ins_mul_mod(px,l,l);
        ins_add(px,px,l);
        if(curve_a2) begin
            //ins_add(px,px)//TODO: ins_inc
            `assert(1,0);
        end
        ins_mul_mod(l,l,px);
        ins_add(py,py,l);
        ins_add(py,py,px);
    end*/
end
endtask

task test_point_double_core;
    input [WORD_WIDTH-1:0] x;
    input [WORD_WIDTH-1:0] y;
    input [WORD_WIDTH-1:0] expected_x;
    input [WORD_WIDTH-1:0] expected_y;
begin
    write_dmem_word(0,x);
    write_dmem_word(1,y);
    ins_mem_ld(0,0);
    ins_mem_ld(1,1);
    point_double();
    check_gpr(0,expected_x);
    check_gpr(1,expected_y);
end
endtask

task test_point_double;
begin
    init_curve_k233();
    test_point_double_core(
        256'h17232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126,
        256'h1db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3,
        256'h1a96a52534c02824c92539163f2ed13243feb57b45adbe4cf7ec61957f6,
        256'h1f9d11ccd5ff37c021bb64dff8df25af3ebc5c3f9bfc5cb17b2203703a8
    );
    $display("test_point_double pass");
end
endtask

initial begin
clk = 0;
forever #5 clk = ~clk;
end

initial begin
$dumpfile ("tb_data_path.fst");
$dumpvars(0,tb_data_path);

test_mem_ops();
test_alu_add();
test_alu_shl(5);
test_mul_plain(2,0,1);
test_mul_mod(2,0,1);
test_div_mod(2,0,1);
test_point_double();
$display("sim pass");
$finish();
end


endmodule
