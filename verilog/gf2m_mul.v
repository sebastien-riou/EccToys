`timescale 1ns / 1ps
`define WORD_WIDTH 256

//implemented:
//op   __-----__
//run  ___----__
//done ______-__
//out  xxxxxxxvv


//other considered option
//op   __-----__
//start__-______
//run  ___----__
//early_____-___
//done ______-__
//out  xxxxxxxvv




module gf2m_mul #(
    parameter MUL_STEP = 8
    )(
    input wire clk,
    input wire reset,
    input wire mod_mul,
    input wire plain_mul,
    input wire stos,
    input wire stox,
    input wire dbus_sel,
    input wire clear,
    input wire [`WORD_WIDTH-1:0] irreducible_poly_msb,
    input wire [`WORD_WIDTH-1:0] irreducible_poly,
    input wire [`WORD_WIDTH-1:0] sbus,
    output reg [`WORD_WIDTH-1:0] dbus,
    output reg run,
    output reg done
    );

reg mul;
always @* mul = mod_mul|plain_mul;
reg run_done;
always @* done = stos|stox|run_done;

reg [10:0] cnt;

reg [`WORD_WIDTH-1:0] s;//scanned operand
reg [`WORD_WIDTH-1:0] x;//xored operand
reg [2*`WORD_WIDTH-1:0] out;

always @* dbus = dbus_sel ? out[`WORD_WIDTH+:`WORD_WIDTH] : out[0+:`WORD_WIDTH];

reg [2*`WORD_WIDTH-1:0] next_out;
reg [`WORD_WIDTH-1:0] next_s;
reg [`WORD_WIDTH-1:0] next_x;

function  [3*`WORD_WIDTH-1:0] mul_step;
    input [`WORD_WIDTH-1:0] s;
    input [`WORD_WIDTH-1:0] x;
    input [2*`WORD_WIDTH-1:0] out;
    reg [`WORD_WIDTH-1:0] next_s;
    reg [2*`WORD_WIDTH-1:0] next_out;
    reg [`WORD_WIDTH-1:0] tmp_s;
    reg [2*`WORD_WIDTH-1:0] tmp_out;
    integer i;
    begin
        next_s = s;
        next_out = out;
        for(i=0;i<MUL_STEP;i=i+1) begin: LOOP
            next_s = {next_s[`WORD_WIDTH-2:0],next_s[`WORD_WIDTH-1]};
            next_out = {next_out[2*`WORD_WIDTH-2:0],1'b0};
            if(next_s[0]) begin
                next_out = next_out ^ x;
            end
            if(mod_mul) begin
                if(next_out & irreducible_poly_msb) begin
                    next_out = next_out ^ irreducible_poly;
                end
            end
        end
        mul_step = {next_s,next_out};
    end
 endfunction

always @* begin
    next_x = x;
    if(stox) next_x = sbus;
end

always @(posedge clk,posedge reset) begin
    if(reset) begin
        s <= 0;
        out <= 0;
    end else begin
        if(mul) begin
            {s,out} <= mul_step(s,x,out);
        end else begin
            if(clear) out <= 0;
            if(stos)  s <= sbus;
        end
    end
end

always @(posedge clk,posedge reset) begin
    if(reset) x <= 0;
    else begin
        x <= next_x;
    end
end


always @(posedge clk,posedge reset) begin
    if(reset) begin
        run <= 0;
        run_done <= 0;
        cnt <= 0;
    end else begin
        if(mul) begin
            run <= 1'b1;
            cnt <= cnt + MUL_STEP;
        end
        if(cnt==(`WORD_WIDTH-MUL_STEP)) begin
            run_done <= 1;
            run <= 1'b0;
            cnt <= 0;
        end
        if(run_done) begin
            run_done <= 1'b0;
        end
    end
end



endmodule
